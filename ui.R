require(shiny)
require(shinythemes)
require(shinyBS)
require(shinydashboard)
require(plotly)
require(d3heatmap)
require(DT)
require(pcaExplorer)
require(manhattanly)
require(shinyFiles)

tweak_group_in <-
  tags$head(tags$style(HTML("
                            .multicol .shiny-options-group{
                            -webkit-column-count: 3; /* Chrome, Safari, Opera */
                            -moz-column-count: 3;    /* Firefox */
                            column-count: 3;
                            -moz-column-fill: balanced;
                            -column-fill: balanced;
                            }
                            .checkbox{
                            margin-top: 0px !important;
                            -webkit-margin-after: 0px !important;
                            }
                            .checkbox-inline {margin: 0 !important;}
                            ")))

landing_pan <- tabPanel("Main",
                        sidebarLayout(
                          sidebarPanel(width = 3,
                                       #upload stuff here
                                       # actionButton("load_user", label = "Load User Dataset", class = "success", width = "250px"),
                                       shinyDirButton("load_user", label = "Load User Dataset", title = NULL, buttonType = "success"),
                                       br(),br(),
                                       actionButton("load_x", label = "Load Example Dataset", width = "200px"),
                                       br(),br(),
                                       radioButtons("corr", "Select Transformation", choices = c("log2", "VST"), selected = "log2"),
                                       bookmarkButton(id = "bookmark_test")
                          ),
                          mainPanel(
                            h2("Welcome to RNAseeker"),
                            h4("Please upload your data using the side bar and begin exploring it"),
                            h5("You can perform"),
                            tags$li("Quality control"),
                            tags$li("Principal component analysis"),
                            tags$li("Exploritory data analysis"),
                            tags$li("Differential expression analysis"),
                            tags$li("Gene ontology enrichment"),
                            br(),
                            textOutput("notes"),
                            br(),
                            # sample_select
                            uiOutput("sel_samp_out"),
                            br(),
                            uiOutput("sel"),
                            DT::dataTableOutput('tbl'),
                            plotlyOutput("mil_reads")
                          )
                        )
)


QC_pan <- tabPanel("QC",
                   mainPanel(
                     #this opend the multiQC file
                     uiOutput("QC")
                   )
)

QC_rrna <- tabPanel("rRNA contamination",
                    mainPanel(
                      h2("rRNA contamination accross sampels"),
                      plotlyOutput("rRNA_cont")
                    )
)


QC_noramlized_hist <- tabPanel("Density Plots",
                               mainPanel(
                                 # h2("Raw Counts"),
                                 # plotlyOutput("dens_raw", width = "100%"),
                                 h2("Density Plot"),
                                 p("This filter will be applied to the datset in the other tabs"),
                                 fluidRow(column(6,uiOutput("gene_count")),
                                          column(6,uiOutput("gene_count_sample"))),
                                 plotlyOutput("dens_log", width = "100%")
                                 # h2("Variance Transformed Counts")
                               )
)




sampdist_pan <- tabPanel("Sample to sample distances and PCA",
                         # sidebarLayout(
                         # sidebarPanel(width = 3,
                         # sliderInput("pcs", label = h4("Number of PCs to retain"),min = 1, max = 10,value = 2)
                         # ),

                         mainPanel(width = 12,
                                   h2("Sample to Sample distances"),
                                   fluidRow(column(3,uiOutput("geneslide")),
                                            column(3,uiOutput("geneslide_box"))
                                   ),
                                   d3heatmapOutput("heatmapsampledist"),

                                   br(), br(),

                                   p("Calculating correlation matricies can be time consuming, click the update button to process and update plot"),

                                   fluidRow(
                                     column(2, uiOutput("col_h")),
                                     column(2, selectInput("pair_corr", "Correlation Method", choices = c("pearson", "kendall", "spearman"), selected = 1)),
                                     column(2, selectInput("k_num", "Number of clusters", choices = c(1:10), selected = 1)),
                                     column(2, selectInput("dist_met", "dist method", choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = 1)),
                                     column(2, selectInput("hclust", "hclust method", choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"), selected = 1)),
                                     column(2, actionButton("corr_up", "Update Corr", class= "btn-primary"))
                                   ),

                                   plotlyOutput("heatmapsampledist_corr"),
                                   h2("PCA plot"),
                                   # sliderInput("pcs", label = h4("Number of PCs to retain"),min = 1, max = 10,value = 2),
                                   fluidRow(
                                     column(4, selectInput("PC1", "X-axis PC", choices = c(1:10), selected = 1)),
                                     column(4, selectInput("PC2", "Y-axis PC", choices = c(1:10), selected = 2)),
                                     column(2, uiOutput("col_pca"))
                                   ),
                                   tabBox(
                                     id = "tabset_pca", #height = "250px",
                                     width = 12,

                                     tabPanel("PCA static with elipses",
                                              plotOutput("pca")
                                     ),

                                     tabPanel("PCA interactive",
                                              plotlyOutput("pca_plotly")
                                     ),

                                     tabPanel("Scree plot",
                                              plotOutput("scree")
                                     ),
                                     tabPanel("3D PCA plot",
                                              p("add selector here"),
                                              plotlyOutput("pca_3d")
                                     )
                                   ),
                                   # actionButton("pca_exp", label = "Launch interactive PCA explorer", class = "success"),
                                   br(), br(),
                                   box(width = 12, height = "100%",
                                       bsCollapse(id = "collapseTest", #open = "Further interactive PCA explorer - click to expand",
                                                  bsCollapsePanel("Further interactive PCA explorer - click to expand", #style = "bigbox",
                                                                  # if (exists("dds_HTSeq"))
                                                                  # pcaExplorer(dds = dds_HTSeq(), rlt = assay_in())
                                                                  pcaExplorer()
                                                  ))),
                                   br(), br(), br(),br(), br(), br(),br(), br(), br()
                         )
)


sat_pan <- tabPanel("Saturation plots",
                    mainPanel(
                      uiOutput("slider_sat_ui"),
                      br(),
                      tabBox(
                        id = "tabset3", #height = "250px",
                        width = 12,
                        tabPanel("Overall",
                                 plotOutput("sat_1")
                        ),
                        tabPanel("Protien coding",
                                 plotOutput("sat_2")
                        )
                      )
                    )
)


sens_pan <- tabPanel("Sensitivity plots",
                     mainPanel(
                       uiOutput("slider_sens_ui"),
                       br(),
                       tabBox(
                         id = "tabset2", #height = "250px",
                         width = 12,
                         tabPanel("Overall",
                                  plotOutput("sens_1")
                         ),
                         tabPanel("Protien coding",
                                  plotOutput("sens_2")
                         )
                       )
                     )
)


comps_pan <- tabPanel("Sample to sample comparison plots",
                      mainPanel(
                        uiOutput("slider_comps_ui"),
                        br(),
                        tabBox(
                          id = "tabset1", #height = "250px",
                          width = 12,
                          tabPanel("Overall",
                                   plotOutput("bio_plot_1")
                          ),
                          tabPanel("Protien coding", plotOutput("bio_plot_2"))
                        )
                      )
)










advanced_options_DE <- box(width = 10, height = "100%",
                           bsCollapse(id = "collapseDEadv",
                                      bsCollapsePanel("More differential expression options", #style = "bigbox",
                                                      fluidRow(
                                                        column(2, sliderInput("lcf", "LCF Threshold", 0, 5, value = 0, step = 0.5)  ),
                                                        column(2, textInput("cook", "Cooks Cutoff", value = "Auto" )  ),
                                                        column(2, textInput("de_beta", "Beta Prior", value = "Auto" )  ),
                                                        column(2, selectInput("p_adjM", "p Adjust Method",
                                                                              choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), selected = ("BH"))),
                                                        bsTooltip("cook", "theshold on Cook's distance, such that if one or more samples for a row have a distance higher, the p-value for the row is set to NA",
                                                                  "right", options = list(container = "body")),
                                                        bsTooltip("p_adjM", "The method to use for adjusting p-values",
                                                                  "right", options = list(container = "body")),
                                                        bsTooltip("de_beta", "whether or not to put a zero-mean normal prior on the non-intercept coefficients See nbinomWaldTest for description of the calculation of the beta prior. By default, the beta prior is used only for the Wald test, but can also be specified for the likelihood ratio test.",
                                                                  "right", options = list(container = "body")),
                                                        bsTooltip("lcf", "Specifies a log2 fold change threshold",
                                                                  "right", options = list(container = "body"))
                                                      ),
                                                      fluidRow(
                                                        # theta
                                                        # the quantiles at which to assess the number of rejections from independent filtering
                                                        column(2, selectInput("alt", "Alt Hypothesis",
                                                                              choices = c("Auto",  "greaterAbs", "lessAbs", "greater", "less"), selected = ("Auto"))),
                                                        column(2, selectInput("de_filt", "Filter Type",
                                                                              choices = c("Auto", "parametric", "local", "mean"), selected = ("Auto"))),
                                                        column(2, selectInput("de_test", "Test to apply",
                                                                              choices = c("Auto", "Wald", "LRT"), selected = ("Auto")) ),

                                                        bsTooltip("de_test", "will use either Wald significance tests, or the likelihood ratio test on the difference in deviance between a full and reduced model formula",
                                                                  "right", options = list(container = "body")),
                                                        bsTooltip("alt", "Those values of log2 fold change which the user is interested in finding. If the log2 fold change specified by name or by contrast is written as beta, then the possible values for altHypothesis represent the following alternate hypotheses: NULL: automatically determine the correct test, greaterAbs: |beta| > lfcThreshold, and p-values are two-tailed lessAbs: |beta| < lfcThreshold, NOTE: this requires that betaPrior=FALSE has been specified in the previous DESeq call. p-values are the maximum of the upper and lower tests. greater: beta > lfcThreshold less: beta < -lfcThreshold",
                                                                  "right", options = list(container = "body")),

                                                        bsTooltip("de_filt", "parametric - fit a dispersion-mean relation of the form: dispersion = asymptDisp + extraPois / mean via a robust gamma-family GLM. The coefficients asymptDisp and extraPois are given in the attribute coefficients of the dispersionFunction of the object. local - use the locfit package to fit a local regression of log dispersions over log base mean (normal scale means and dispersions are input and output for dispersionFunction). The points are weighted by normalized mean count in the local regression. mean - use the mean of gene-wise dispersion estimates.",
                                                                  "right", options = list(container = "body"))
                                                      )
                                      )))
#add to server code
#results
# lfcThreshold = input$lcf, , pAdjustMethod = input$p_adjM, cooksCutoff = input$cook,
#dds
# fitType=input$de_filt, betaPrior = input$betaprior




subplots_volc <- tabBox(id = "volc_sub_box", #height = "250px",
                        width = 12,
                        tabPanel("Box Plot",
                                 plotlyOutput("boxplot_sel")
                                 # uiOutput("boxplot_sel")
                        ),
                        tabPanel("Heatmap",
                                 plotlyOutput("heat_sel")
                        ),
                        tabPanel("test",
                                 verbatimTextOutput("brush")
                        )
)


de_pan <- tabPanel("Differential Expression",
                   mainPanel(
                     fluidRow(
                       column(3, sliderInput("alpha_in", "Set alpha", value = 0.05, step = 0.01, max = 0.2, min = 0.001)),
                       column(2, textInput("p_in",label="Set p filter", value = 0.1)),
                       column(4, uiOutput("meta_deUI"))
                     ),
                     # p("More differential expression options"),
                     h3("Add a button to save the results as a pdf and also save the r object and to load that r object!"),
                     advanced_options_DE,
                     fluidRow(
                       column(3, div(uiOutput("de_grp_1_UI"))),
                       column(3, div(uiOutput("de_grp_2_UI"))),
                       column(3, tableOutput("sample_table"))
                     ),
                     actionButton("DEcalc", "Differential Expression Calculation"),

                     bsTooltip("DEcalc", "Click here to perform the DE analysis, this may take a while..",
                               "right", options = list(container = "body")),

                     br(), br(),

                     tabBox(
                       id = "tabset_de", #height = "250px",
                       width = 12,
                       tabPanel("Table",
                                dataTableOutput("DE_stat")
                       ),
                       tabPanel("Volcano",
                                plotlyOutput("DE_volc")
                       ),
                       tabPanel("QQ",
                                plotlyOutput("DE_qq")
                       ),
                       tabPanel("MA",
                                plotlyOutput("DE_ma")
                       ),
                       tabPanel("Scatter",
                                # p("Plotly not working for this at the mo"),
                                plotlyOutput("DE_scat")
                       )
                     ),
                     br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                     textOutput("genes_sel"),
                     br(), br(),
                     subplots_volc
                     # verbatimTextOutput("console_DE")
                   )
)



ge_go <- tabPanel("GO enrichment",
                  mainPanel(
                    h4("Gene set enrichment"),
                    p("still under construction"),
                    textInput("p_go", "Enter adj.p cutoff for GO results (alpha)", value = 0.1),
                    selectInput("method_go", "GO enrichment method",
                                choices = c("anova", "randomForest"), selected = ("anova")),
                    bsTooltip("method_go", "The statistical framework to score genes and gene ontologies."),
                    actionButton("Calculate", "Calculate"),

                    verbatimTextOutput("Go_progress"),
                    p("Genes regarded as significant using the alpha set in DE panel will be included in the analysis"),
                    plotOutput("pwd_plot"),
                    verbatimTextOutput("enriched"),
                    tableOutput("res_tab")
                  ))

function(request) {
  fluidPage(
    theme = shinythemes::shinytheme("slate"),
    # includeCSS("/Users/jdlim/Downloads/x_lbd_free_v1.3.1/assets/css/bootstrap.min.css"),
    tweak_group_in,
    # theme = "www/slate.edited.css",
    # tags$style(type = "text/css", "#bigbox {height: 300px !important;}"),
    navbarPage(title=div(img(src = "test_logo_white.png", height="30", width="30"), "RNAseq analysis"),
               windowTitle = "RNAseeker - SB",

               landing_pan,

               navbarMenu("Quality Control",
                          QC_rrna,
                          QC_pan,
                          QC_noramlized_hist),

               navbarMenu("Visual Exploration",
                          # QC_noramlized_hist,
                          sat_pan,
                          sens_pan,
                          comps_pan,
                          sampdist_pan),

               navbarMenu("Statistical Exploration",
                          de_pan,
                          ge_go
               )

    )
  )
  # )
}