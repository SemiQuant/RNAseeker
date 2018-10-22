# using install.packages(c("ggplot2", "shiny", "plotly"))
#other devtools::install_github("rstudio/shiny"); devtools::install_github("hadley/ggplot2"); devtools::install_github("ropensci/plotly")

# rsconnect::deployApp("/Users/jdlim/Bioinformatics/RNAseeker", launch.browser = F, account = "sciencebuff")

#to fix
# updating of corr plots after sample selection

#########
##TO DO##
#########
# make seletion on volcan plot show other plots - have hover and select as radiobutton options
#
# check scatter plot of counts is correct
#
# add GO analysis
#
# add 'googleAuthR' button and to load data
# add bookmark button
#
# add annotation info to plotly and tables
# add annotation lookup
#
# make plots look good
#
# add download as HTML widget for plots
#
#################################################################################################
source("qqly.R")
require(shiny)
require(ggplot2)
require(plotly)
require(DESeq2)
require(NOISeq)
require(adegenet)
require(DT)
require(ggfortify)
require(heatmaply)
require(BiocParallel)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(goseq)
require(GO.db)
require(GOexpress)


#custom theme for ggplot
gg_back <- theme(
  panel.background = element_rect(fill = "#272b30"),
  plot.background = element_rect(fill = "#272b30"),
  legend.background = element_rect(fill = "#272b30"),
  panel.grid.major = element_line(colour = "black"),
  axis.title.x = element_text(colour = "light grey"),
  axis.title.y = element_text(colour ="light grey"),
  legend.title = element_text(colour = "white"),
  plot.title = element_text(colour = "light grey"),
  legend.text = element_text(colour = "light grey"),
  axis.text = element_text(colour ="light grey")
)

testing <<- "T"
threads <<- 1
register(MulticoreParam(threads))


shinyServer(function(input, output, session) {
  #############
  ##load data##
  #############
  data_load <- reactiveValues()

  observeEvent(input$load_x,{
    data_load$rRNA_file_pat <- "example_data/rRNA_check/"
    data_load$vst_path <- "example_data/vst_in"
    data_load$rlt_path <- "example_data/rlt_in"
    data_load$dds_path <-"example_data/dds_HTSeq_in"
    data_load$bio_path <- "bio.txt"
    data_load$noi_sat_path <- "example_data/noi_dat_saturation"
    data_load$noi_coutns_path <- "example_data/noi_dat_countsbio"
    data_load$noi_bio_path <- "example_data/noi_dat_bio_detect"
    data_load$multi_qc_path <- "/Users/jdlim/Library/Mobile Documents/com~apple~CloudDocs/Bioinformatics/RNAseeker/example_data/multiqc_report.html"
  })


  observeEvent(input$load_user,{
    data_load$rRNA_file_pat <- data_load$vst_path <- data_load$rlt_path <- data_load$dds_path <- data_load$bio_path <- data_load$noi_sat_path <- data_load$noi_coutns_path <- data_load$noi_bio_path <- data_load$multi_qc_path <- NULL

    data_load$rRNA_file_pat <- paste0(input$load_user, "/rRNA_check/")
    data_load$vst_path <- paste0(input$load_user, "/vst_in")
    data_load$rlt_path <- paste0(input$load_user, "/rlt_in")
    data_load$dds_path <-paste0(input$load_user, "/dds_HTSeq_in")
    data_load$bio_path <- paste0(input$load_user, "bio.txt")
    data_load$noi_sat_path <- paste0(input$load_user, "/noi_dat_saturation")
    data_load$noi_coutns_path <- paste0(input$load_user, "/noi_dat_countsbio")
    data_load$noi_bio_path <- paste0(input$load_user, "/noi_dat_bio_detect")
    data_load$multi_qc_path <- paste0(input$load_user, "/multiqc_report.html")

  #   data_load$snp_matrix <- as.character(input$meta_dat$datapath[1])
  #   data_load$genp_path <- as.character(input$snp_mat_IN$datapath[1])
  })

  rlt_in <- reactive({
    load(data_load$rlt_path)
    if (testing == "T")
      rlt_in[1:500]
    else
      rlt_in
  })

    dds_HTSeq_in <- reactive({
      load(data_load$dds_path)
      if (testing == "T")
        dds_HTSeq_in[1:500]
      else
        dds_HTSeq_in
    })

      vst_in <- reactive({
        load(data_load$vst_path)
        if (testing == "T")
          vst_in[1:500]
        else
          vst_in
      })


      #   # bio_detect <- eventReactive(input$load_x,{
      #   #   load("example_data/bio_detect")
      #   #   names(bio_detect@dat$biotables) <- gsub(".Homo_sapiens.HTSeq.counts", "", names(bio_detect@dat$biotables))
      #   #   bio_detect
      #   # }, ignoreInit = TRUE)
      #   #
  # rlt_in <- eventReactive(input$load_x, input$load_u,{
  #   print(data_load$rlt_path)
  #   load(data_load$rlt_path)
  #   if (testing == "T")
  #     rlt_in[1:500]
  #   else
  #     rlt_in
  # })


  noi_dat <- reactive({
    bio <- read.table(data_load$bio_path)
    readData(assay(dds_HTSeq_in()), colData(dds_HTSeq_in()), biotype = bio)
  })

  noi_dat_saturation <- reactive({
    # noi_dat_saturation() <- dat(noi_dat(), k = 0, ndepth = 10, type = "saturation")
    load(data_load$noi_sat_path)
    noi_dat_saturation
  })

  noi_dat_countsbio <- reactive({
    # noi_dat_countsbio() <- dat(noi_dat(), factor = NULL, type = "countsbio")
    load(data_load$noi_coutns_path)
    noi_dat_countsbio
  })

  noi_dat_bio_detect <- reactive({
    # noi_dat_bio_detect() <- dat(noi_dat(), k = 0, type = "biodetection", factor = NULL)
    load(data_load$noi_bio_path)
    noi_dat_bio_detect
  })

  # observeEvent(input$load_x, input$load_u,{
  #   data_load$multi_qc_path <<- "/Users/jdlim/Library/Mobile Documents/com~apple~CloudDocs/Bioinformatics/RNAseeker/example_data/multiqc_report.html"
  # })


  #load data end#
  ###############

  #create the main tabke with samples (this is where samples can be selected)
  sa_tab <- eventReactive(input$load_x,{
    s_tab <- data.frame(colData(dds_HTSeq_in()))
    s_tab$sample <- row.names(s_tab)
    conts <- data.frame(colSums(assay(dds_HTSeq_in())))
    s_tab <- merge(s_tab, conts, by = "row.names")
    genes <- data.frame(colSums(assay(dds_HTSeq_in())>1))
    sa_tab <- merge(s_tab, genes, by.x = "sample", by.y = "row.names")[,2:6]
    colnames(sa_tab) <- c("Sample", "Group", "Rep", "Read Count", "Features Detected")
    sa_tab
  })


  observeEvent(c(input$load_x, input$load_user),{

    # render radio buttons based on the groups in the dds object to be used for selection in table
    output$sel <- renderUI({
      all_samples <- colnames(dds_HTSeq_in())
      all_groups <- as.character(unique(colData(dds_HTSeq_in())[,1]))
      tags$div(align = 'left',
               class = 'multicol',
               # checkboxGroupInput("sel_samp", "Select samples to keep by name",
               #                    inline = T,
               #                    choiceNames = as.list(all_samples),
               #                    choiceValues = as.list(all_samples),
               #                    selected = as.list(all_samples)
               # ),
               checkboxGroupInput("sel_samp_g", "Select samples to keep by group (not working yet)",
                                  inline = T,
                                  # choiceNames = as.list(all_groups),
                                  # choiceValues = as.list(all_groups),
                                  choices = as.list(all_groups),
                                  selected = as.list(all_groups)
               )
      )
    })

    output$sel_samp_out <- renderUI(actionButton("sel_samp", h4("Update Sample Selection")))

    #data table of samples, metadata from dds and reads counts - see sa_tab function
    output$tbl <- DT::renderDataTable({
      datatable(sa_tab(),
                selection = list(target = 'row', selected = 1:nrow(sa_tab()))
      ) %>% formatStyle(
        c(colnames(sa_tab())),
        backgroundColor = "red",
        color = "black"
      )
    },  server = TRUE)

    output$gene_count <- renderUI({
      sliderInput("min_count", "Keep features with more than this many normalized counts", min = 0, max = 100, value = 5, step = 1)
    })

    output$gene_count_sample <- renderUI({
      sliderInput("min_count_sample", "In this many samples", min = 1, max = length(input$tbl_rows_selected), value = length(input$tbl_rows_selected), step = 1)
    })

    output$mil_reads <- renderPlotly({
      counts_total <- counts(dds_HTSeq_in())
      counts_total <- colSums(counts_total)
      counts_total <- data.frame(counts_total)
      counts_total$Sample <- rownames(counts_total)

      c <- ggplot(counts_total, aes(Sample,counts_total)) +
        geom_bar(stat = "identity", aes(fill = counts_total)) + #colour = counts_total,
        # scale_colour_gradient(low = "blue", high = "green") +
        scale_fill_gradientn(colours = heat.colors(10), name = "Count\nIntensity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        gg_back +
        ggtitle("Counts of reads per sample") +
        ylab("Count") + xlab(NULL)
      ggplotly(c) #%>% layout(autosize = FALSE, width = "400px")
    })


  }, ignoreInit = T)


  selected_samples <- reactive({
    input$sel_samp
    sa_tab()[input$tbl_rows_selected,]$Sample
  })

  dds_HTSeq <- reactive({
    input$sel_samp
    dds_HTSeq <- dds_HTSeq_in()[,selected_samples()]
    dds_HTSeq <- estimateSizeFactors(dds_HTSeq)
    idx <- try(rowSums( counts(dds_HTSeq, normalized=TRUE) >= input$min_count ) >= input$min_count_sample)

    if (class(idx) == "try-error")
      dds_HTSeq
    else
      dds_HTSeq[idx, ]
  })

  assay_in <- reactive({
    input$sel_samp
    if (input$corr == "log2")
      assay_in <- rlt_in()
    else
      assay_in <- vst_in()

    assay_in[rownames(dds_HTSeq()), colnames(assay_in)%in%selected_samples()]
  })


  expgroups <- reactive({
    input$load_x
    input$sel_samp
    expgroups <- as.data.frame(colData(dds_HTSeq()))
    rownames(expgroups) <- colnames(dds_HTSeq())
    expgroups
  })

  output$notes <- renderText(readLines("example_data/user_notes.txt"))

  output$dens_log <- renderPlotly({
    # input$sel_samp
    # de_counts <- stack(de_counts)
    # de_counts <- de_counts[de_counts$values >0, ]
    # log2(values+1)

    #or raw counts
    # counts(dds_HTseq)
    req(c(input$min_count, input$min_count_sample), cancelOutput = T)

    #this is for the selected assay (vst or rlog)
    # de_counts <- assay(assay_in())
    # de_counts <- stack(de_counts)
    # de_counts <- data.frame(de_counts[c(2,4)])
    # colnames(de_counts) <- c("Sample", "value")
    # g <- ggplot(de_counts, aes(value, colour = Sample, fill = Sample )) +
    #   geom_density(alpha = 0.05) + theme(legend.position="none") +
    #   xlab(input$corr)


    #this is for raw, normalized counts
    de_counts <- counts(dds_HTSeq(), normalized = T)
    de_counts <- stack(de_counts)
    de_counts <- data.frame(de_counts[c(2,4)])
    colnames(de_counts) <- c("Sample", "value")
    g <- ggplot(de_counts, aes(value, colour = Sample, fill = Sample )) +
      geom_density(alpha = 0.05) + theme(legend.position="none") +
      xlab("Log10 of normalized counts") + scale_x_log10() +
      ylab("Density")

    g <- g + gg_back

    ggplotly(g)

  })



  #plot rRNA contamination
  observeEvent(c(input$load_x, input$load_user, input$sel_samp),{
    if (file.exists(data_load$rRNA_file_pat)){
      #rRNA
      files_rrna <- list.files(data_load$rRNA_file_pat, full.names = T)
      rrna <- NULL
      for (file in files_rrna){
        tmp.r <- read.table(file)
        rrna <- rbind(rrna, tmp.r[c(6,9)])
      }
      rrna$V9 <- gsub("%","",rrna$V9)
      rrna$V9 <- as.numeric(rrna$V9)
      rrna$V6 <- gsub("_R1_.*", "", rrna$V6)
      colnames(rrna) <- c("Sample", "Percentage_rRNA")

      rrna <- rrna[rrna$Sample%in%selected_samples(),]

      r_plot <- ggplot(data = rrna, aes(x = Sample, y = Percentage_rRNA)) +
        geom_histogram(stat = "identity") + #, aes(fill = Sample)
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylab("Percentage rRNA reads") + xlab(NULL)


      r_plot <- r_plot + gg_back

      #rRNA end
      output$rRNA_cont <- renderPlotly({ggplotly(r_plot)})
    }else{
      output$rRNA_cont <- renderPlotly({
        ggplotly(
          ggplot(data.frame()) +
            annotate("text", x=8, y=13000, label= "No Data", size = 20, color = "red") +
            gg_back
        )
      })
    }
  }, ignoreInit = T)
  #plot rRNA contamination end

  #not working??
  observeEvent(c(input$load_x, input$load_user),
               output$QC <- renderUI(browseURL(data_load$multi_qc_path, encodeIfNeeded = T))
  )

  #sample dists
  output$geneslide <- renderUI({
    sliderInput("g_slide", "Select number of top genes to retain", min = 2, max = length(rownames(dds_HTSeq())), value =  length(rownames(dds_HTSeq())))
  })

  output$geneslide_box <- renderUI({
    textInput("g_box", "", value =  length(rownames(dds_HTSeq())))
  })

  observeEvent(input$g_box,{
    updateSliderInput(session, "g_slide", value = input$g_box)
  })

  assay_red <- reactive({
    input$g_slide
    mads <- apply(assay(assay_in()), 1, mad)

    #or can use the row variation like this
    # library(genefilter)
    # ntop <- input$g_slide
    # rv <- rowVars(assay(assay.tmp))
    # select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    #
    # assay(assay.tmp)[select, ]

    assay_in()[order(mads, decreasing=T)[1:input$g_slide], ]
  })




  output$heatmapsampledist <- renderD3heatmap({
    # if (!is.null(input$color_by)) {
    # assay_red.tmp <- assay_red()
    # expgroups() <- as.data.frame(colData(assay_red.tmp)[, "Group"])
    # rownames(expgroups()) <- colnames(assay_red.tmp)
    # colnames(expgroups()) <- "Group"
    d3heatmap(as.matrix(dist(t(assay(assay_red())))),
              annotation_col = expgroups())
  })


  output$col_h <- renderUI({
    selectInput("h_col", "Color by metadata", choices = names(colData(dds_HTSeq_in())), selected = 1)
  })

  corr <- reactive({
    input$corr_up
    isolate(cor(assay(assay_red()), method = input$pair_corr))
  })

  observe({
    output$heatmapsampledist_corr <- renderPlotly({
      # if(!is.null(corr()))
      heatmaply(corr(),
                row_side_colors = fac2col(expgroups()[,colnames(expgroups())==input$h_col]),
                col_side_colors = fac2col(expgroups()[,colnames(expgroups())==input$h_col]),
                k_row = as.numeric(input$k_num), k_col = as.numeric(input$k_num),
                hclust_method = input$hclust, dist_method = input$dist_met,
                scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limits = c(-1, 1)),
                heatmap_layers = gg_back,
                side_color_layers = gg_back
      )

    })
  })

  output$col_pca <- renderUI({
    selectInput("pc_col", "Color by metadata", choices = names(colData(dds_HTSeq_in())), selected = 1)
  })

  output$pca <- renderPlot({
    # expgroups() <- as.data.frame(colData(dds_HTSeq_in()))
    # rownames(expgroups()) <- colnames(dds_HTSeq_in())
    # colnames(expgroups()) <- "Group"

    pc_as <- t(assay(assay_red()))
    # pc_as <- t(assay(rlt_in))
    pca <- prcomp(pc_as, scale. = F, center = T)
    eig <- pca$sdev^2

    output$pca_plotly <- renderPlotly(
      autoplot(pca, data = expgroups(), x = input$PC1, y = input$PC2,
               colour = input$pc_col, frame.colour = input$pc_col) + #frame = TRUE, frame.type = 'norm'
        gg_back
    )

    output$scree <- renderPlot({
      screeplot(pca,
                col = "#d3d3d3", col.lab = "#d3d3d3", col.main = "white", col.sub = "white", col.axis = "#d3d3d3", col.axis = "#d3d3d3", main = "")
      }, bg = "#272b30")

    output$pca_3d <- renderPlotly({
      # d3.tmp <- data.frame(pca$rotation)
      # d3.tmp <- merge(d3.tmp, expgroups()[c(input$pc_col)], by = "row.names")
      # rownames(d3.tmp) <- d3.tmp$Row.names; d3.tmp <- d3.tmp[,-1]
      #
      # plot_ly(d3.tmp, x = ~PC1, y = ~PC2, z = ~PC3, color = input$pc_col) %>%
      #   add_markers() %>%
      #   layout(scene = list(xaxis = list(title = 'PC1'),
      #                       yaxis = list(title = 'PC2'),
      #                       zaxis = list(title = 'PC3')))


      ###3d scatter param
      axis <- list(
        xaxis = list(title = 'PC1'),
        yaxis = list(title = 'PC2'),
        zaxis = list(title = 'PC3'),
        color =toRGB("#d3d3d3"),
        linecolor = toRGB("#d3d3d3"),
        gridcolor = toRGB("#d3d3d3")
      )

      scene = list(
        xaxis = axis,
        yaxis = axis,
        zaxis = axis)
      ####

      if (all(rownames(expgroups()[c(input$pc_col)]) == rownames(data.frame(pca$x))))
      plot_ly(data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color = expgroups()$Group,
              colors = unique(fac2col(grps$group)), alpha = 0.8) %>%
        add_markers() %>%
        layout(scene = scene,
               paper_bgcolor = "#272b30",
               plot_bgcolor = "#272b30",
               legend = list(
                 font = list(
                   color = "#d3d3d3")
               )
        )
        # plot_ly(data.frame(pca$rotation), x = ~PC1, y = ~PC2, z = ~PC3, color = unlist(expgroups()[c(input$pc_col)])) %>%
        # add_markers() %>%
        # layout(scene = list(xaxis = list(title = 'PC1'),
        #                     yaxis = list(title = 'PC2'),
        #                     zaxis = list(title = 'PC3')))
      else
        ggplotly(
          ggplot(data.frame()) +
            annotate("text", x=8, y=13000, label= "No Data", size = 20, color = "red") +
            gg_back
        )
    })

    # pca <- dudi.pca(pc_as, center=TRUE, scale=FALSE, scannf = F, nf = input$pcs)
    #cool pca plots
    # par(bg = "grey")
    # s.class(pca$l1, expgroups()$Group, col = unique(fac2col(expgroups()$Group)))


    # s.class(pca$x[,c(input$PC1, input$PC2)], expgroups()[,input$pc_col], col = unique(fac2col(expgroups()[,input$pc_col])))
    s.class(data.frame(pca$x), xax = as.numeric(input$PC1), yax = as.numeric(input$PC2),
            expgroups()[,colnames(expgroups())==input$pc_col], col = unique(fac2col(expgroups()[,colnames(expgroups())==input$pc_col])))
    add.scatter.eig(eig, 3, 1, c(as.numeric(input$PC1), as.numeric(input$PC2)), posi = "bottomright")
                    # bg.col = "#272b30")





  }, bg = "#272b30")

  #sat, sens and sample plots
  observeEvent(c(input$load_x, input$load_user),{

    output$slider_sat_ui <- renderUI({

      #change this to the sample names
      selectInput("slider_sat", "Select samples to compare",
                  choices = selected_samples(),
                  multiple = T,
                  width = "100%",
                  selected = selected_samples()[1:3])

      # sliderInput(inputId = "slider_sat", "Select sample range to explore", min = 1, max = length(names(noi_dat_saturation()@dat$depth)), value = c(1, 10), step = 1)
      #length(input$tbl_rows_selected)
      # input$tbl_rows_selected
    })

    output$sat_1<- renderPlot({
      par(col.lab="white", col="grey", col.axis = "white", col.main = "white", col.sub = "white")
      # explo.plot(noi_dat_saturation(), toplot = 1, samples = c(input$slider_sat))
      explo.plot(noi_dat_saturation(), toplot = 1, samples = c(1:length(names(noi_dat_saturation()@dat$depth)))[names(noi_dat_saturation()@dat$depth)%in%c(input$slider_sat)],
                 col.ticks = "white")
    }, bg = "#272b30")

    output$sat_2 <- renderPlot({
      par(col.lab="white", col="grey", col.axis = "white", col.main = "white", col.sub = "white")
      explo.plot(noi_dat_saturation(), toplot = "protein_coding", samples = c(1:length(names(noi_dat_saturation()@dat$depth)))[names(noi_dat_saturation()@dat$depth)%in%c(input$slider_sat)],
                 col.ticks = "white")
    }, bg = "#272b30")

    output$slider_sens_ui <- renderUI({
      selectInput("slider_sens", "Select samples to compare",
                  choices = selected_samples(),
                  multiple = T,
                  width = "100%",
                  selected = selected_samples()[1:3])
    })

    output$sens_1 <- renderPlot({
      par(col.lab="white", col="grey", col.axis = "white", col.main = "white", col.sub = "white")
      explo.plot(noi_dat_countsbio(), toplot = 1, plottype = "barplot",
                 samples = c(1:length(names(noi_dat_countsbio()@dat$summary$global)[!(names(noi_dat_countsbio()@dat$summary$global)%in%c("global", "total"))]))[names(noi_dat_countsbio()@dat$summary$global)[!(names(noi_dat_countsbio()@dat$summary$global)%in%c("global", "total"))]%in%c(input$slider_sens)],
                 col.ticks = "white")
    }, bg = "#272b30")

    output$sens_2 <- renderPlot({
      par(col.lab="white", col="grey", col.axis = "white", col.main = "white", col.sub = "white")
      explo.plot(noi_dat_countsbio(), toplot = "protein_coding", plottype = "boxplot",
                 samples = c(1:length(names(noi_dat_countsbio()@dat$summary$global)[!(names(noi_dat_countsbio()@dat$summary$global)%in%c("global", "total"))]))[names(noi_dat_countsbio()@dat$summary$global)[!(names(noi_dat_countsbio()@dat$summary$global)%in%c("global", "total"))]%in%c(input$slider_sens)],
                 col.ticks = "white")
    }, bg = "#272b30")

    output$slider_comps_ui <- renderUI({
      fluidRow(
        column(4,

               selectInput("slider_comps1", "Select samples to compare",
                           choices = selected_samples(),
                           multiple = F,
                           width = "100%",
                           selected = selected_samples()[1])),
        column(4,
               selectInput("slider_comps2", "Select samples to compare",
                           choices = selected_samples(),
                           multiple = F,
                           width = "100%",
                           selected = selected_samples()[2]))
      )
    })

    output$bio_plot_1 <- renderPlot({
      par(col.lab="white", col="grey", col.axis = "white", col.main = "white", col.sub = "white")
      par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
      explo.plot(noi_dat_bio_detect(), plottype = "persample", samples = c(1:length(names(noi_dat_bio_detect()@dat$biotables)))[names(noi_dat_bio_detect()@dat$biotables)%in%c(input$slider_comps1, input$slider_comps2)],
                 col.ticks = "white")
    }, bg = "#272b30")

    output$bio_plot_2 <- renderPlot({
      par(col.lab="white", col="grey", col.axis = "white", col.main = "white", col.sub = "white")
      par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
      explo.plot(noi_dat_bio_detect(), toplot = "protein_coding", plottype = "comparison", samples = c(1:length(names(noi_dat_bio_detect()@dat$biotables)))[names(noi_dat_bio_detect()@dat$biotables)%in%c(input$slider_comps1, input$slider_comps2)],
                 col.ticks = "white")
    }, bg = "#272b30")

    output$bio_plot_2_text <- renderText({
      print("Show test for difference btw samples here")
      #   # par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
      #   explo.plot(noi_dat_bio_detect(), toplot = "protein_coding", plottype = "comparison", samples = c(1:length(names(noi_dat_bio_detect()@dat$biotables)))[names(noi_dat_bio_detect()@dat$biotables)%in%c(input$slider_comps1, input$slider_comps2)])
    })

  }, ignoreInit = T)


  #statistical tab
  # withProgress(message = 'Performing Differential Expression Calculations',
  #              detail = 'This may take a while...', value = 0, {
  #                # incProgress(1/2)
  #              })


  #need selector with
  output$meta_deUI <- renderUI(selectInput("meta_de", "Select Variable Group to compare",
                                           choices = colnames(expgroups()),
                                           selected = 1))

  # then radio buttons with
  output$de_grp_1_UI <- renderUI({
    radioButtons("de_grp_1", "Group 1",
                 choices = unique(unlist(expgroups()[,input$meta_de])),
                 selected = unique(unlist(expgroups()[,input$meta_de]))[1])
  })

  # then radio buttons with
  output$de_grp_2_UI <- renderUI({
    radioButtons("de_grp_2", "Group 2", choices = unique(unlist(expgroups()[,input$meta_de])),
                 selected = unique(unlist(expgroups()[,input$meta_de]))[2])
  })

  output$sample_table <- renderTable(table(expgroups()[, input$meta_de]))

  observeEvent(input$DEcalc, {
    req(input$de_grp_1, input$de_grp_2)

    isolate({
      if(input$de_grp_1 !=input$de_grp_2){

        # a character vector with exactly three elements:
        # the name of a factor in the design formula,
        # the name of the numerator level for the fold change,
        # and the name of the denominator level for the fold change (simplest case)

        de_options <- NULL
        if (input$de_test != "Auto")
          de_options <- paste0(de_options, 'test="', input$de_test, '",')
        if (input$de_filt != "Auto")
          de_options <- paste0(de_options, 'fitType="', input$de_filt, '",')
        if (input$de_beta != "Auto")
          de_options <- paste0(de_options, 'betaPrior="', input$de_beta, '",')


        # if (testing == "dfd"){
        #   # save(dds, "example_data/testing/dds")
        #   load("example_data/testing/dds")
        #   dds<-dds
        #   # save(rld.test, "example_data/testing/rld.test")
        #   load("example_data/testing/rld.test")
        #   rld.test<-rld.test
        #   # save(counts.test, "example_data/testing/counts.test")
        #   load("example_data/testing/counts.test")
        #   counts.test<-counts.test
        #   # save(counts.all, "example_data/testing/counts.all")
        #   load("example_data/testing/counts.all")
        #   counts.all<-counts.all
        #   # save(res, "example_data/testing/res")
        #   load("example_data/testing/res")
        #   res<-res
        #   # save(dds.test, "example_data/testing/dds.test")
        #   load("example_data/testing/dds.test")
        #   dds.test<-dds.test
        #
        # }else{

        dds <- reactive({
          dds.test <- dds_HTSeq()
          dds.test <- dds.test[,colData(dds.test)[,colnames(colData(dds.test))==input$meta_de]%in%c(input$de_grp_1, input$de_grp_2)]
          colData(dds.test)[,colnames(colData(dds.test))==input$meta_de] <- as.factor(as.character(colData(dds.test)[,colnames(colData(dds.test))==input$meta_de]))
          dds.test <- estimateSizeFactors(dds.test)

          if (is.null(de_options))
            dds <- DESeq(dds.test, parallel = T)
          else
            dds <- eval(parse(text =
                                cat(
                                  paste("DESeq(dds.test, parallel = T,",
                                        de_options,
                                        ")"))))
          dds
        })
        # rld.test <- rlog(dds, blind = FALSE) #for hetamps and such
        counts.test <- log2(counts(dds(), normalized=T))
        counts.all <- log2(counts(dds_HTSeq(), normalized=T))

        res <- reactive({
          res_opts <- NULL
          if (input$alt != "Auto")
            res_opts <- paste0(res_opts, 'altHypothesis="', input$alt, '",')
          if (input$cook != "Auto")
            res_opts <- paste0(res_opts, 'cooksCutoff="', input$cook, '",')

          if (is.null(res_opts))
            res <- results(dds(),
                           alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2),
                           pAdjustMethod = input$p_adjM, lfcThreshold = input$lcf, parallel = T)
          else
            res <- eval(parse(text =
                                cat(
                                  paste("results(dds, parallel = T, alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2), pAdjustMethod = input$p_adjM,lfcThreshold = input$lcf,",
                                        res_opts,
                                        ")"))))

          #Annotate
          res$symbol <- mapIds(org.Hs.eg.db,
                               keys=row.names(res),
                               column="SYMBOL",
                               keytype="ENSEMBL",
                               multiVals="first")
          res$entrez <- mapIds(org.Hs.eg.db,
                               keys=row.names(res),
                               column="ENTREZID",
                               keytype="ENSEMBL",
                               multiVals="first")

          res

        })

        # if (input$de_beta == "Auto")
        #   dds <- DESeq(dds_HTSeq(), test = input$de_test, fitType = input$de_filt)
        # else
        #   dds <- DESeq(dds_HTSeq(), test = input$de_test, fitType = input$de_filt, betaPrior = input$de_beta)
        #
        # if (input$cook == "Auto")
        #   res <- results(dds,
        #                   alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2),
        #                   altHypothesis = input$alt, pAdjustMethod = input$p_adjM, lfcThreshold = input$lcf,
        #                   test = input$de_test)
        # else
        #   res <- results(dds,
        #                   alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2),
        #                   altHypothesis = input$alt, pAdjustMethod = input$p_adjM, lfcThreshold = input$lcf,
        #                   test = input$de_test,
        #                   cooksCutoff = input$cook)

        # parallel = T
        # library(BiocParallel)
        # register(MulticoreParam(10))
        # DESeq(dds, parallel=10)
        # results(DESeq_dds, parallel = T)

        #DE outputs
        output$DE_stat <- renderDataTable({
          # req(input$de_grp_1, input$de_grp_2)
          res <- res()
          res <- res[order(res$padj, decreasing = F),]
          # tmp <- data.frame(res[1:sum(res$padj < input$p_in, na.rm=TRUE), ])
          tmp <- data.frame(res[!(is.na(res$pvalue)),])
          tmp[,1:6] <- round(tmp[,1:6], 3)
          tmp$Gene <- row.names(tmp)
          tmp <- tmp[c(9,1,2,3,4,5,6,7,8)]
          datatable(tmp, rownames = F) %>% formatStyle(
            'padj',
            backgroundColor = styleInterval(input$p_in, c('#98FB98', '#D3D3D3')),
            color = "black"
          )  %>% formatStyle(
            c(colnames(tmp)),
            color = "black"
          )
        })

        res_ma <- data.frame(res())
        res_ma$Fold_change <- ifelse(res_ma$log2FoldChange <0, "down", "up")
        res_ma$Fold_change <- ifelse(is.na(res_ma$Fold_change), "none", res_ma$Fold_change)
        res_ma$Significat_at_p_adjusted <- ifelse(res_ma$padj <= input$p_in, "yes", "no")


        output$DE_ma <- renderPlotly({
          # plotMA(res, ylim=c(-5,5))
          ggplotly(
            ggplot(res_ma, aes(x=baseMean, y=log2FoldChange)) +
              geom_point(aes( colour = Significat_at_p_adjusted, shape=Fold_change )) + scale_x_continuous(trans='log10') +
              geom_hline(yintercept = 0, col = "red") + xlab("Mean of Normalized Counts") + ylab("Log Fold Change") +
              ylim(c(-5, 5)) +
              geom_smooth(se = F) +
              gg_back
          )

        })

        ### volc
        # input=NULL
        # input$p_in=0.1
        volc <- data.frame(res())
        volc$P <- volc$pvalue
        volc$log_pvalue <- -log10(volc$pvalue)
        volc$sig_p.adj <- ifelse(volc$padj <= input$p_in, "yes", "no")
        volc$Gene <- paste0(rownames(volc), " (", volc$symbol, ")")
        volc$key_test <- seq(1:nrow(volc))

        output$DE_volc <- renderPlotly({
          vol <- ggplot(data=volc,
                        aes(x = log2FoldChange, y = log_pvalue, colour = sig_p.adj, text = Gene, key = key_test)) +
            geom_point(alpha=0.4, size=1.75) +
            # xlim(c(-6, 6)) +
            xlab("log2FoldChange") + ylab("-log10 p-value (unadjusted)") +
            theme_bw() +
            theme(legend.position="right")+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            gg_back

          ggplotly(vol, tooltip = c("log2FoldChange", "log_pvalue", "Gene"), source = "volc_sor") %>% layout(dragmode = "select")
          # %>% layout(dragmode = "select")

          #try changing to plot_ly as it seems to be much less sluggish
          # plot_ly(volc, x = ~log2FoldChange, y = ~ log_pvalue,
          #         source = "volc_sor")
        })

        output$DE_qq <- renderPlotly({
          qqly(na.omit(volc))
        })


        output$DE_scat <- renderPlotly({

          #test
          # load("example_data/delete/dds_2")
          # dds
          # rld.test <- rlog(dds, blind = FALSE)
          # counts.test <- counts(dds, normalized=T)
          # res <- results(dds)
          #
          #
          # res[1,]
          #
          # assay(rld.test)[1,]
          # counts.test[1,]
          #
          #
          #
          # scat.tmp <- assay(rld.test)
          # scat.tmp.1 <- data.frame(scat.tmp[1,colnames(scat.tmp)%in%rownames(colData(rld.test)[colData(rld.test)$Group == "LTBI",])])
          # scat.tmp.2 <- data.frame(scat.tmp[1,colnames(scat.tmp)%in%rownames(colData(rld.test)[colData(rld.test)$Group == "prev",])])
          # scat.tmp <- NULL
          # cond1 <- colMeans(scat.tmp.1, na.rm = T)
          # cond2 <- colMeans(scat.tmp.2, na.rm = T)
          #
          #
          # res[1,]
          # cond2-cond1
          #
          #
          #
          #
          #
          # scat.tmp <- counts.test
          # scat.tmp.1 <- data.frame(scat.tmp[1,colnames(scat.tmp)%in%rownames(colData(rld.test)[colData(rld.test)$Group == "LTBI",])])
          # scat.tmp.2 <- data.frame(scat.tmp[1,colnames(scat.tmp)%in%rownames(colData(rld.test)[colData(rld.test)$Group == "prev",])])
          # scat.tmp <- NULL
          # cond1 <- colMeans(scat.tmp.1, na.rm = T)
          # cond2 <- colMeans(scat.tmp.2, na.rm = T)
          #
          #
          # res[1,]
          # cond2-cond1
          # log2(cond2)-log2(cond1)
          #

          #test end

          #should recalculate this using the subsetted dds
          # scat.tmp <- assay(rld.test)
          expgroups.tmp <- expgroups()[c(colnames(expgroups())==input$meta_de)]
          colnames(expgroups.tmp) <- "Group"
          expgroups.tmp$sample <- rownames(expgroups.tmp)

          scat.tmp <- counts.test
          scat.tmp.1 <- data.frame(scat.tmp[,colnames(scat.tmp)%in%
                                              expgroups.tmp[expgroups.tmp$Group==input$de_grp_1, ]$sample
                                            ])
          scat.tmp.2 <- data.frame(scat.tmp[,colnames(scat.tmp)%in%
                                              expgroups.tmp[expgroups.tmp$Group==input$de_grp_2, ]$sample
                                            ])
          scat.tmp <- NULL
          scat.tmp.1$cond1 <- rowMeans(scat.tmp.1, na.rm = T)
          scat.tmp.2$cond2 <- rowMeans(scat.tmp.2, na.rm = T)
          scat <- data.frame("cond1" = scat.tmp.1$cond1,
                             "cond2" = scat.tmp.2$cond2,
                             "Gene" = rownames(scat.tmp.1),
                             "Significat_at_p_adjusted" = res_ma$Significat_at_p_adjusted,
                             "Adj.p" = res_ma$padj)

          p.scat <- ggplot(scat, aes(x=cond1, y=cond2)) +
            geom_point(aes( colour = Significat_at_p_adjusted)) +
            xlab(input$de_grp_1) + ylab(input$de_grp_2) + gg_back

          #jerryrig to get extra in the tiptool..
          p.scat <- p.scat + geom_point(aes(Adj.p), alpha = 0) +
            geom_text(aes(label=Gene), hjust=0, vjust=0, alpha = 0)

          ggplotly(p.scat, tooltip = c("Gene", "cond1", "cond2")) #, "Adj.p"
        })


        output$brush <- renderPrint({
          d <- event_data("plotly_selected", source = "volc_sor", session)
          if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
        })


        observe({
          #setup tmp data for subplots
          eve_data <- event_data("plotly_selected", source = "volc_sor", session)
          volc_key <- as.numeric(eve_data$key)
          if (!(is.null(eve_data))){
            # if (length(volc_key) == 1)
            #   count.tmp <- t(data.frame(counts.test[c(volc_key),]))
            # else
            #   count.tmp <- counts.test[c(volc_key),]
            count.tmp <- counts.test[c(volc_key),, drop=F]

            count.tmp <- data.frame(t(count.tmp))
            count.tmp$mean_count <- rowMeans(count.tmp)
            ano <- colData(dds())[,input$meta_de]
            # table(rownames(colData(dds.test))==rownames(count.tmp))
            #all samples
            # if (length(volc_key) == 1)
            #   counts.all.tmp <- t(data.frame(counts.all[c(volc_key),]))
            # else
            #   counts.all.tmp <- counts.all[c(volc_key),]
            counts.all.tmp <- counts.all[c(volc_key),, drop=F]

            counts.all.tmp <- data.frame(t(counts.all.tmp))
            counts.all.tmp$mean_count <- rowMeans(counts.all.tmp)
            counts.all.tmp$Sample <- rownames(counts.all.tmp)
            counts.all.tmp$Genes <- paste(colnames(counts.all.tmp)[1:(length(colnames(counts.all.tmp))-2)], sep="", collapse=",")

            ano.all <- colData(dds_HTSeq())[,input$meta_de]

            output$genes_sel <- renderText(c("Seleceted Genes", colnames(count.tmp)[colnames(count.tmp) != "mean_count"] ))

            output$boxplot_sel <- renderPlotly({
              p1 <- ggplotly(
                ggplot(count.tmp, aes(x = ano, y = mean_count, color = ano)) +
                  scale_y_log10() + geom_boxplot(outlier.alpha = 0) + geom_jitter(cex = 3) +
                  xlab(input$meta_de) + ylab("Mean count accross selected genes")
              )
              # no need to have its own data here, just use the highlight() from plotly...?

              p2 <- ggplotly(
                ggplot(counts.all.tmp, aes(x = ano.all, y = mean_count, color = ano.all)) +
                  scale_y_log10() + geom_boxplot(outlier.alpha = 0) + geom_jitter(cex = 3) +
                  xlab(input$meta_de) + ylab("Mean count accross selected genes - all samples in initial table")
              )
              # p2 <- p2 + tooltip = c("Sample", "mean_count", "Genes")

              subplot(p1, p2, shareX = F, shareY = T)

              # }else{
              #   ggplotly(
              #     ggplot(data.frame()) +
              #       annotate("text", x=8, y=13000, label= "Please select some points", size = 20, color = "orange") +
              #       gg_back
              #   )
              # }
              # }
            })

            output$heat_sel <- renderPlotly({
              if (length(volc_key) > 1){
                # count.heat.tmp <- assay(rld.test)[c(volc_key),]
                count.heat.tmp <- counts.test[c(volc_key),]
                # h1 <-
                heatmaply(count.heat.tmp,
                          # col_side_colors = fac2col( expgroups()[c(colnames(expgroups())==input$meta_de)] ),
                          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", mid = "black", high = "green", midpoint = 0, na.value = "white"),
                          xlab = "Mean count accross selected genes",
                          heatmap_layers = gg_back,
                          side_color_layers = gg_back
                )
                #           # k_row = as.numeric(input$k_num), k_col = as.numeric(input$k_num),
                #           # hclust_method = input$hclust, dist_method = input$dist_met,
                #
                # plot_ly(z = count.heat.tmp, type = "heatmap",
                #         col_side_colors = fac2col(colData(dds.test)[,input$meta_de]))

                # tooltip = c("Sample", "mean_count", "Genes")

                # subplot(h1, h2, shareX = F, shareY = T)
                # h1
              }
            })

          }
        })

        ###

        # #this must take whats selected in manhattan plot
        # output$DE_gene <- renderPlotly({
        #   req(input$gene_in)
        #   d <- plotCounts(dds, gene=input$gene_in, intgroup=input$meta_de,
        #                   returnData=TRUE)
        #   p <- ggplot(d, aes(x=condition.Group, y=count)) +
        #     geom_point(position=position_jitter(w=0.1, h=0))
        #   ggplotly()
        # })


        ############
        ##go panel##
        ############
        # require(goseq)
        observeEvent(input$Calculate, {
          load("hg19.ensGene.LENGTH")
          res.go <- res()
          res.go <- res.go[!is.na(res.go$pvalue),]
          sig_genes <- ifelse(res.go$padj<=input$p_in, 1, 0)
          names(sig_genes) <- rownames(res.go)
          # need to download the annotations
          pwf <- nullp(sig_genes, "hg19", "ensGene", plot.fit = F)

          output$pwd_plot <- renderPlot({
            plotPWF(pwf)
          })

          output$enriched <- renderText({
            GO.wall <- goseq(pwf, "hg19", "ensGene")
            enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<input$p_go]
            for(go in enriched.GO){
              print(GOTERM[[go]])
              cat("--------------------------------------\n")
            }
          })


          counts.tmp <- log2(counts(dds(), normalized=T))
          counts.phen <- annotatedDataFrameFrom(as.matrix(colData(dds())), byrow = T)

          counts.go <- ExpressionSet(counts.tmp, counts.phen)
          pData(counts.go)$Group <- colData(dds())$Group

          #download annotations
          # listDatasets(useMart('ensembl'))
          # listMarts()
          # ensembl89 = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
          #
          # allgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'),
          #                          mart=ensembl89)
          # colnames(allgenes.Ensembl)[1] = 'gene_id'
          # allGO.Ensembl = getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'),
          #                       mart=ensembl89)
          # allGO.Ensembl = allGO.Ensembl[allGO.Ensembl$go_id != '',]
          #
          # GOgenes.Ensembl = getBM( attributes=c('ensembl_gene_id', 'go_id'),
          #                          mart=ensembl89)
          #
          # colnames(GOgenes.Ensembl)[1] = 'gene_id'
          #
          # GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
          # GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]
          #
          # save(GOgenes.Ensembl, file = "GOgenes.Ensembl")
          # save(allGO.Ensembl, file = "allGO.Ensembl")
          # save(allgenes.Ensembl, file = "allgenes.Ensembl")
          load("GOgenes.Ensembl")
          load("allGO.Ensembl")
          load("allgenes.Ensembl")

          set.seed(1987)
          GO_results <- GO_analyse(eSet = counts.go, f = "Group",
                                   GO_genes=GOgenes.Ensembl, all_GO=allGO.Ensembl, all_genes=allgenes.Ensembl,
                                   method = input$method_go) #rf
          #rf not working for some reason

          # output$Go_progress <- renderText({

          #add progress meter
          GO_results.pVal <- pValue_GO(result=GO_results)
          # })
          output$res_tab <- renderTable({
            data.frame(GO_results.pVal$GO)
          })

          # need to make assignment for input$GO_selected, liked to table?
          output$GO_heat <- renderPlot({
            heatmap_GO(input$GO_selected, GO_results, counts.go)
          })

        })

      }
    })
  }, ignoreInit = T) #end DE action button


  observeEvent(input$bookmark_test, {
    session$doBookmark()
  })



  #
  # observeEvent(input$test_save,{
  #   save(session, file = "delete/test")
  # })
  #
  # observeEvent(input$test_load,{
  #   session = load("delete/test")
  # })
  #

  # observeEvent(input$test_load,{
  #
  #   if(!file.exists("delete/inputs.RDS")) {return(NULL)}
  #
  #   savedInputs <- readRDS("delete/inputs.RDS")
  #
  #   inputIDs <- names(savedInputs)
  #   inputvalues <- unlist(savedInputs)
  #   for (i in 1:length(savedInputs)) {
  #     session$sendInputMessage(inputIDs[i],  list(value=inputvalues[[i]]) )
  #   }
  # })
  #
  # observeEvent(input$test_save,{
  #   saveRDS( reactiveValuesToList(input) , file = "delete/inputs.RDS")
  # })







  # #for DE calcs
  # save.image(file = "test.R")
  # dds(), counts.test, counts.all, res()
  # onBookmark(function(state) {
  #   state$values$DE <- vals$sum
  # })
  #


})

enableBookmarking(store = "server") #or "url"






