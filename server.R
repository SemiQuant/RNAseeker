# using install.packages(c("ggplot2", "shiny", "plotly"))
#other devtools::install_github("rstudio/shiny"); devtools::install_github("hadley/ggplot2"); devtools::install_github("ropensci/plotly")

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

shinyServer(function(input, output, session) {

  testing = "T"
  threads = 1
  register(MulticoreParam(threads))

  #load data
  observeEvent(input$load_x,{
    #set these by upload, aslo rRNA_file_pat its just commented for testing
    # print("testing")
    multi_qc_path <<- "/Users/jdlim/Library/Mobile Documents/com~apple~CloudDocs/Bioinformatics/RNAseeker/example_data/multiqc_report.html"
    rRNA_file_pat <<- "example_data/rRNA_check/"
    bio <<- read.table("bio.txt")


    # dds_HTSeq
    load("example_data/dds_HTSeq_in")
    dds_HTSeq_in <<- dds_HTSeq_in #make it accessable to UI
    # rlt
    load("example_data/rlt_in")
    rlt_in <<- rlt_in
    # vst
    load("example_data/vst_in")
    vst_in <<- vst_in


    load("example_data/bio_detect")
    names(bio_detect@dat$biotables) <- gsub(".Homo_sapiens.HTSeq.counts","",names(bio_detect@dat$biotables))
    bio_detect <<- bio_detect

    # assay_in <- vst_in
    # assay_red <- assay_in
    #calc both here and then load later, will be quicker then redoing each time like now
    # mads <- apply(assay(assay_in), 1, mad)

    #initilize - do same for user data
    # vstd <- vst_in
    # rlt <- rlt_in
    # dds_HTSeq <- dds_HTSeq_in
    selected_samples <<- rownames(colData(dds_HTSeq_in))
    expgroups <<- as.data.frame(colData(dds_HTSeq_in))
    rownames(expgroups) <<- colnames(dds_HTSeq_in)


    noi_dat <<- readData(assay(dds_HTSeq_in), colData(dds_HTSeq_in), biotype = bio)
    # noi_dat_saturation <<- dat(noi_dat, k = 0, ndepth = 10, type = "saturation")
    load("example_data/noi_dat_saturation")
    noi_dat_saturation <<- noi_dat_saturation

    # noi_dat_countsbio <<- dat(noi_dat, factor = NULL, type = "countsbio")
    load("example_data/noi_dat_countsbio")
    noi_dat_countsbio <<- noi_dat_countsbio

    # noi_dat_bio_detect <<- dat(noi_dat, k = 0, type = "biodetection", factor = NULL)
    load("example_data/noi_dat_bio_detect")
    noi_dat_bio_detect <<- noi_dat_bio_detect





    if (testing == "T"){
      rlt_in <<- rlt_in[1:50]
      vst_in <<- vst_in[1:50]
      dds_HTSeq_in <<- dds_HTSeq_in[1:50]
    }



  })


  observeEvent(c(input$load_x, input$load_user),{

    output$sel <- renderUI({
      all_samples <<- colnames(dds_HTSeq_in)
      all_groups <- as.character(unique(colData(dds_HTSeq_in)[,1]))
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
                                  choiceNames = as.list(all_groups),
                                  choiceValues = as.list(all_groups),
                                  selected = as.list(all_groups)
               )
      )
    })

    output$sel_samp_out <- renderUI({actionButton("sel_samp", h4("Update Sample Selection"))})

    output$tbl <- DT::renderDataTable({
      s_tab <- data.frame(colData(dds_HTSeq_in))
      s_tab$sample <- row.names(s_tab)
      conts <- data.frame(colSums(assay(dds_HTSeq_in)))
      s_tab <- merge(s_tab, conts, by = "row.names")
      genes <- data.frame(colSums(assay(dds_HTSeq_in)>1))
      sa_tab <<- merge(s_tab, genes, by.x = "sample", by.y = "row.names")[,2:6]
      colnames(sa_tab) <<- c("Sample", "Group", "Rep", "Read Count", "Features Detected")
      datatable(sa_tab,
                selection = list(target = 'row', selected = 1:nrow(sa_tab)) )

    }, server = TRUE)

    output$gene_count <- renderUI({
      sliderInput("min_count", "Keep features with more than this many normalized counts", min = 0, max = 100, value = 5, step = 1)
    })

    output$gene_count_sample <- renderUI({
      sliderInput("min_count_sample", "In this many samples", min = 1, max = length(input$tbl_rows_selected), value = length(input$tbl_rows_selected), step = 1)
    })


    output$mil_reads <- renderPlotly({
      counts_total <- counts(dds_HTSeq_in)
      counts_total <- colSums(counts_total)
      counts_total <- data.frame(counts_total)
      counts_total$Sample <- rownames(counts_total)

      c <- ggplot(counts_total, aes(Sample,counts_total)) +
        geom_bar(stat = "identity", aes(fill = counts_total)) + #colour = counts_total,
        # scale_colour_gradient(low = "blue", high = "green")+
        scale_fill_gradientn(colours = heat.colors(10)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        gg_back +
        ggtitle("Counts of reads per sample")
      ggplotly(c) #%>% layout(autosize = FALSE, width = "400px")
    })


  }, ignoreInit = T)


  observeEvent(input$sel_samp, {
    selected_samples <<- sa_tab[input$tbl_rows_selected,]$Sample
  }, ignoreInit = T)

  dds_HTSeq <<- reactive({
    input$sel_samp
    dds_HTSeq <- dds_HTSeq_in[,selected_samples]
    dds_HTSeq <- estimateSizeFactors(dds_HTSeq_in)
    idx <- try(rowSums( counts(dds_HTSeq, normalized=TRUE) >= input$min_count ) >= input$min_count_sample)

    if (class(idx) == "try-error")
      dds_HTSeq
    else
      dds_HTSeq[idx, ]
  })

  assay_in <<- reactive({
    if (input$corr == "log2")
      assay_in <- rlt_in
    else
      assay_in <- vst_in

    assay_in[rownames(dds_HTSeq()), colnames(assay_in)%in%selected_samples]
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
      xlab("Log10 of normalized counts") + scale_x_log10()

    g <- g + gg_back

    ggplotly(g)

  })





  #plot rRNA contamination
  observeEvent(c(input$load_x, input$load_user, input$sel_samp),{
    if (exists("rRNA_file_pat")){
      #rRNA
      files_rrna <- list.files(rRNA_file_pat, full.names = T)
      rrna <- NULL
      for (file in files_rrna){
        tmp.r <- read.table(file)
        rrna <- rbind(rrna, tmp.r[c(6,9)])
      }
      rrna$V9 <- gsub("%","",rrna$V9)
      rrna$V9 <- as.numeric(rrna$V9)
      rrna$V6 <- gsub("_R1_.*", "", rrna$V6)
      colnames(rrna) <- c("Sample", "Percentage_rRNA")

      rrna <- rrna[rrna$Sample%in%selected_samples,]

      r_plot <- ggplot(data = rrna, aes(x = Sample, y = Percentage_rRNA)) +
        geom_histogram(stat = "identity") + #, aes(fill = Sample)
        theme(axis.text.x = element_text(angle = 90, hjust = 1))


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
  })
  #plot rRNA contamination end

  #not working??
  observeEvent(c(input$load_x, input$load_user),
               output$QC <- renderUI(browseURL(multi_qc_path, encodeIfNeeded = T))
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
  #
  #   observeEvent(input$g_slide,{
  #     updateTextInput(session, "g_box", value = input$g_slide)
  #   })



  assay_red <- reactive({
    input$g_slide
    assay.tmp <- assay_in()
    mads <- apply(assay(assay.tmp), 1, mad)

    #or can use the row variation like this
    # library(genefilter)
    # ntop <- input$g_slide
    # rv <- rowVars(assay(assay.tmp))
    # select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    #
    # assay(assay.tmp)[select, ]

    assay.tmp[order(mads, decreasing=T)[1:input$g_slide], ]
  })




  output$heatmapsampledist <- renderD3heatmap({
    # if (!is.null(input$color_by)) {
    # assay_red.tmp <- assay_red()
    # expgroups <- as.data.frame(colData(assay_red.tmp)[, "Group"])
    # rownames(expgroups) <- colnames(assay_red.tmp)
    # colnames(expgroups) <- "Group"
    d3heatmap(as.matrix(dist(t(assay(assay_red())))),
              annotation_col = expgroups)
  })


  output$col_h <- renderUI({
    selectInput("h_col", "Color by metadata", choices = names(colData(dds_HTSeq_in)), selected = 1)
  })

  corr <- reactive({
    input$corr_up
    isolate(cor(assay(assay_red()), method = input$pair_corr))
  })

  observe({
    output$heatmapsampledist_corr <- renderPlotly({
      # if(!is.null(corr()))
      heatmaply(corr(),
                row_side_colors = fac2col(expgroups[,colnames(expgroups)==input$h_col]),
                col_side_colors = fac2col(expgroups[,colnames(expgroups)==input$h_col]),
                k_row = as.numeric(input$k_num), k_col = as.numeric(input$k_num),
                hclust_method = input$hclust, dist_method = input$dist_met,
                scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limits = c(-1, 1))
      )

    })
  })




  output$col_pca <- renderUI({
    selectInput("pc_col", "Color by metadata", choices = names(colData(dds_HTSeq_in)), selected = 1)
  })


  output$pca <- renderPlot({
    # expgroups <- as.data.frame(colData(dds_HTSeq_in))
    # rownames(expgroups) <- colnames(dds_HTSeq_in)
    # colnames(expgroups) <- "Group"

    pc_as <- t(assay(assay_red()))
    # pc_as <- t(assay(rlt_in))

    pca <- prcomp(pc_as, scale. = F, center = T)
    eig <- pca$sdev^2


    output$pca_plotly <- renderPlotly(
      autoplot(pca, data = expgroups, x = input$PC1, y = input$PC2,
               colour = input$pc_col, frame.colour = input$pc_col) + #frame = TRUE, frame.type = 'norm'
        gg_back
    )

    output$scree <- renderPlot(screeplot(pca))


    output$pca_3d <- renderPlotly({
      # d3.tmp <- data.frame(pca$rotation)
      # d3.tmp <- merge(d3.tmp, expgroups[c(input$pc_col)], by = "row.names")
      # rownames(d3.tmp) <- d3.tmp$Row.names; d3.tmp <- d3.tmp[,-1]
      #
      # plot_ly(d3.tmp, x = ~PC1, y = ~PC2, z = ~PC3, color = input$pc_col) %>%
      #   add_markers() %>%
      #   layout(scene = list(xaxis = list(title = 'PC1'),
      #                       yaxis = list(title = 'PC2'),
      #                       zaxis = list(title = 'PC3')))

      if (all(rownames(expgroups[c(input$pc_col)]) == rownames(data.frame(pca$rotation))))
        plot_ly(data.frame(pca$rotation), x = ~PC1, y = ~PC2, z = ~PC3, color = unlist(expgroups[c(input$pc_col)])) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'PC1'),
                            yaxis = list(title = 'PC2'),
                            zaxis = list(title = 'PC3')))
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
    # s.class(pca$l1, expgroups$Group, col = unique(fac2col(expgroups$Group)))

    # s.class(pca$x[,c(input$PC1, input$PC2)], expgroups[,input$pc_col], col = unique(fac2col(expgroups[,input$pc_col])))
    s.class(data.frame(pca$x), xax = as.numeric(input$PC1), yax = as.numeric(input$PC2),
            expgroups[,colnames(expgroups)==input$pc_col], col = unique(fac2col(expgroups[,colnames(expgroups)==input$pc_col])))
    add.scatter.eig(eig, 3, 1, c(as.numeric(input$PC1), as.numeric(input$PC2)))

  })


  #sat, sens and sample plots
  observeEvent(c(input$load_x, input$load_user),{

    output$slider_sat_ui <- renderUI({

      #change this to the sample names
      selectInput("slider_sat", "Select samples to compare",
                  choices = selected_samples,
                  multiple = T,
                  width = "100%",
                  selected = selected_samples[1:3])


      # sliderInput(inputId = "slider_sat", "Select sample range to explore", min = 1, max = length(names(noi_dat_saturation@dat$depth)), value = c(1, 10), step = 1)
      #length(input$tbl_rows_selected)
      # input$tbl_rows_selected
    })

    output$sat_1<- renderPlot({
      # explo.plot(noi_dat_saturation, toplot = 1, samples = c(input$slider_sat))
      explo.plot(noi_dat_saturation, toplot = 1, samples = c(1:length(names(noi_dat_saturation@dat$depth)))[names(noi_dat_saturation@dat$depth)%in%c(input$slider_sat)])
    }, bg = "#272b30")

    output$sat_2 <- renderPlot({
      explo.plot(noi_dat_saturation, toplot = "protein_coding", samples = c(1:length(names(noi_dat_saturation@dat$depth)))[names(noi_dat_saturation@dat$depth)%in%c(input$slider_sat)])
    }, bg = "#272b30")



    output$slider_sens_ui <- renderUI({
      selectInput("slider_sat", "Select samples to compare",
                  choices = selected_samples,
                  multiple = T,
                  width = "100%",
                  selected = selected_samples[1:3])
    })

    output$sens_1 <- renderPlot({
      explo.plot(noi_dat_countsbio, toplot = 1, plottype = "barplot",
                 samples = c(1:length(names(noi_dat_countsbio@dat$summary$global)[!(names(noi_dat_countsbio@dat$summary$global)%in%c("global", "total"))]))[names(noi_dat_countsbio@dat$summary$global)[!(names(noi_dat_countsbio@dat$summary$global)%in%c("global", "total"))]%in%c(input$slider_sat)])
    }, bg = "#272b30")

    output$sens_2 <- renderPlot({
      explo.plot(noi_dat_countsbio, toplot = "protein_coding", plottype = "boxplot",
                 samples = c(1:length(names(noi_dat_countsbio@dat$summary$global)[!(names(noi_dat_countsbio@dat$summary$global)%in%c("global", "total"))]))[names(noi_dat_countsbio@dat$summary$global)[!(names(noi_dat_countsbio@dat$summary$global)%in%c("global", "total"))]%in%c(input$slider_sat)])
    }, bg = "#272b30")



    output$slider_comps_ui <- renderUI({
      fluidRow(
        column(4,

               selectInput("slider_comps1", "Select samples to compare",
                           choices = selected_samples,
                           multiple = F,
                           width = "100%",
                           selected = selected_samples[1])),
        column(4,
               selectInput("slider_comps2", "Select samples to compare",
                           choices = selected_samples,
                           multiple = F,
                           width = "100%",
                           selected = selected_samples[2]))
      )
    })


    output$bio_plot_1 <- renderPlot({
      par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
      explo.plot(noi_dat_bio_detect, plottype = "persample", samples = c(1:length(names(noi_dat_bio_detect@dat$biotables)))[names(noi_dat_bio_detect@dat$biotables)%in%c(input$slider_comps1, input$slider_comps2)])

    }, bg = "#272b30")

    output$bio_plot_2 <- renderPlot({
      par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
      explo.plot(noi_dat_bio_detect, toplot = "protein_coding", plottype = "comparison", samples = c(1:length(names(noi_dat_bio_detect@dat$biotables)))[names(noi_dat_bio_detect@dat$biotables)%in%c(input$slider_comps1, input$slider_comps2)])
    }, bg = "#272b30")

    output$bio_plot_2_text <- renderText({
      print("Show test for difference btw samples here")
      #   # par(mfrow = c(1,2))  # we need this instruction because two plots (one per sample) will be generated
      #   explo.plot(noi_dat_bio_detect, toplot = "protein_coding", plottype = "comparison", samples = c(1:length(names(noi_dat_bio_detect@dat$biotables)))[names(noi_dat_bio_detect@dat$biotables)%in%c(input$slider_comps1, input$slider_comps2)])
    })

  }, ignoreInit = T)




  #statistical tab
  # withProgress(message = 'Performing Differential Expression Calculations',
  #              detail = 'This may take a while...', value = 0, {
  #                # incProgress(1/2)
  #              })




  #need selector with
  output$meta_deUI <- renderUI(selectInput("meta_de", "Select Variable Group to compare",
                                           choices = colnames(expgroups),
                                           selected = 1))

  # then radio buttons with
  output$de_grp_1_UI <- renderUI({
    radioButtons("de_grp_1", "Group 1",
                 choices = unique(unlist(expgroups[,input$meta_de])),
                 selected = unique(unlist(expgroups[,input$meta_de]))[1])
  })

  # then radio buttons with
  output$de_grp_2_UI <- renderUI({
    radioButtons("de_grp_2", "Group 2", choices = unique(unlist(expgroups[,input$meta_de])),
                 selected = unique(unlist(expgroups[,input$meta_de]))[2])
  })

  output$sample_table <- renderTable(table(expgroups[, input$meta_de]))











  observeEvent(input$DEcalc, {
    req(input$de_grp_1, input$de_grp_2)

    isolate({
      if(input$de_grp_1 !=input$de_grp_2){

        # a character vector with exactly three elements:
        #   the name of a factor in the design formula,
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
        #   dds<<-dds
        #   # save(rld.test, "example_data/testing/rld.test")
        #   load("example_data/testing/rld.test")
        #   rld.test<<-rld.test
        #   # save(counts.test, "example_data/testing/counts.test")
        #   load("example_data/testing/counts.test")
        #   counts.test<<-counts.test
        #   # save(counts.all, "example_data/testing/counts.all")
        #   load("example_data/testing/counts.all")
        #   counts.all<<-counts.all
        #   # save(res, "example_data/testing/res")
        #   load("example_data/testing/res")
        #   res<<-res
        #   # save(dds.test, "example_data/testing/dds.test")
        #   load("example_data/testing/dds.test")
        #   dds.test<<-dds.test
        #
        # }else{

        dds.test <- dds_HTSeq()
        dds.test <- dds.test[,colData(dds.test)[,colnames(colData(dds.test))==input$meta_de]%in%c(input$de_grp_1, input$de_grp_2)]
        colData(dds.test)[,colnames(colData(dds.test))==input$meta_de] <- as.factor(as.character(colData(dds.test)[,colnames(colData(dds.test))==input$meta_de]))
        dds.test <<- estimateSizeFactors(dds.test)

        if (is.null(de_options))
          dds <<- DESeq(dds.test, parallel = T)
        else
          dds <<- eval(parse(text =
                               cat(
                                 paste("DESeq(dds.test, parallel = T,",
                                       de_options,
                                       ")"))))

        # print("dds done")

        rld.test <<- rlog(dds, blind = FALSE) #for hetamps and such
        counts.test <<- counts(dds, normalized=T)
        counts.all <<- counts(dds_HTSeq(), normalized=T)

        res_opts <- NULL
        if (input$alt != "Auto")
          res_opts <- paste0(res_opts, 'altHypothesis="', input$alt, '",')
        if (input$cook != "Auto")
          res_opts <- paste0(res_opts, 'cooksCutoff="', input$cook, '",')

        if (is.null(res_opts))
          res <<- results(dds,
                          alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2),
                          pAdjustMethod = input$p_adjM, lfcThreshold = input$lcf, parallel = T)
        else
          res <<- eval(parse(text =
                               cat(
                                 paste("results(dds, parallel = T, alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2), pAdjustMethod = input$p_adjM,lfcThreshold = input$lcf,",
                                       res_opts,
                                       ")"))))

        # }


        # if (input$de_beta == "Auto")
        #   dds <- DESeq(dds_HTSeq(), test = input$de_test, fitType = input$de_filt)
        # else
        #   dds <- DESeq(dds_HTSeq(), test = input$de_test, fitType = input$de_filt, betaPrior = input$de_beta)
        #
        # if (input$cook == "Auto")
        #   res <<- results(dds,
        #                   alpha=input$alpha_in, contrast = c(input$meta_de, input$de_grp_1, input$de_grp_2),
        #                   altHypothesis = input$alt, pAdjustMethod = input$p_adjM, lfcThreshold = input$lcf,
        #                   test = input$de_test)
        # else
        #   res <<- results(dds,
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
          # res <- res()
          res <- res[order(res$padj, decreasing = F),]
          tmp <- data.frame(res[1:sum(res$padj < input$p_in, na.rm=TRUE), ])
          tmp <- round(tmp,3)
          tmp$Gene <- row.names(tmp)
          tmp[c(7,1,2,3,4,5,6)]
        })

        res_ma <- data.frame(res)
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

        volc <- data.frame(res)
        volc$P <- volc$pvalue
        volc$log_pvalue <- -log10(volc$pvalue)
        volc$sig_p.adj <- ifelse(volc$padj <= input$p_in, "yes", "no")
        volc$Gene <- rownames(volc)
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


        output$brush <- renderPrint({
          d <- event_data("plotly_selected", source = "volc_sor", session)
          if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
        })




        observe({
          #setup tmp data for subplots
          eve_data <- event_data("plotly_selected", source = "volc_sor", session)
          volc_key <- as.numeric(eve_data$key)
          if (!(is.null(eve_data))){
            if (length(volc_key) == 1)
              count.tmp <- t(data.frame(counts.test[c(volc_key),]))
            else
              count.tmp <- counts.test[c(volc_key),]
            count.tmp <- data.frame(t(count.tmp))
            count.tmp$mean_count <- rowMeans(count.tmp)
            ano <- colData(dds.test)[,input$meta_de]
            # table(rownames(colData(dds.test))==rownames(count.tmp))
            #all samples
            if (length(volc_key) == 1)
              counts.all.tmp <- t(data.frame(counts.all[c(volc_key),]))
            else
              counts.all.tmp <- counts.all[c(volc_key),]
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
                count.heat.tmp <- assay(rld.test)[c(volc_key),]
                # h1 <-
                heatmaply(count.heat.tmp,
                          col_side_colors = fac2col(rep(colData(dds.test)[,input$meta_de], 5) ),
                          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", mid = "black", high = "green", midpoint = 0, na.value = "white"),
                          xlab = "Mean count accross selected genes")
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




        # output$DE_scat <- plotlyOutput({
        #
        #   #should recalculate this using the subsetted dds
        #   scat.tmp <- assay(rld.test)
        #   scat.tmp.1 <- data.frame(scat.tmp[,colnames(scat.tmp)%in%rownames(expgroups[expgroups$Group==input$de_grp_1,])])
        #   scat.tmp.2 <- data.frame(scat.tmp[,colnames(scat.tmp)%in%rownames(expgroups[expgroups$Group==input$de_grp_2,])])
        #   scat.tmp <- NULL
        #   scat.tmp.1$cond1 <- rowMeans(scat.tmp.1, na.rm = T)
        #   scat.tmp.2$cond2 <- rowMeans(scat.tmp.2, na.rm = T)
        #   scat <- data.frame("cond1" = scat.tmp.1$cond1,
        #                      "cond2" = scat.tmp.2$cond2,
        #                      row.names = rownames(scat.tmp.1),
        #                      Significat_at_p_adjusted = res_ma$Significat_at_p_adjusted)
        #   ggplotly(
        #     ggplot(scat, aes(x=cond1, y=cond2)) +
        #       geom_point(aes( colour = Significat_at_p_adjusted )) +
        #       xlab(input$de_grp_1) + ylab(input$de_grp_2)
        #   )
        # })


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





      }
    })
  }, ignoreInit = T) #end de action button









})





