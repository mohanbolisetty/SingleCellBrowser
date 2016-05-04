
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)

shinyUI(navbarPage(
  a('SingleCellBiology', href = 'https://www.jax.org/research-and-faculty/tools/scientific-research-services/genome-tech-single-cell-biology/single-cell-biology', target = '_blank',
    style = 'color: black;'),
  windowTitle = 'Single Cell Biology',
  
  tabPanel('overview',
           fluidRow(
             div(h3('Single Cell Live'), align = 'center')
           ),
           fluidRow(
             column(10, offset = 1,
                    p('This app is designed for exploratory data analysis of
                       processed RNA-Seq data of single cell experiments. There are four menu tabs
                      that can be used to choose plots and tables to view.'),
                    
                    p(strong('Information:')),
                    tags$ul(
                      tags$li(strong('v0.0.1'), ': Test'),
                      tags$li(strong('Experiment'), ':T-Cells were sorted, and simulated. Cells were captured on C1'),
                      tags$li(strong('Processing'), ':Alignments GRCh38, ensembl v70. Filtered on Protein coding and lincRNAs')
                    )
                    )),
           br(),
           br(),
           fluidRow(
             column(5,offset=3,
                    sliderInput("cpm_threshold",
                                "Filter Genes < log2CPM:",
                                min=1,
                                max=20,
                                value=2)
                    ),
             column(5, offset=3,
                   sliderInput("cells_threshold",
                               "No. of cells with log2CPM:",
                               min=1,
                               max=20,
                               value=5)
                   )
             )
           ),
  
  navbarMenu('Data QC',
  
          tabPanel('Primary Data QC',
                   fluidRow(
                     div(h4('Data QC - Mapping and Library metrics'),align='center')
                   ),
                   fluidRow(
                     column(3,
                            selectInput('hm_qc', label = 'Plot Type:', choices = c('Total Reads','Mapped Reads','Duplication Rates','Percent Target'), selected = 'Total Reads')
                     )
                    ),
                   fluidRow(plotOutput('qc_plot')),
                   fluidRow(uiOutput("download_hm_plt_button"),'Download Plot')
                   ),
          
          tabPanel('Secondary Data QC',
                   fluidRow(
                     div(h4('Data QC - Genes Detected and Mitochondrial Alignment'),align='center')
                   ),
                   fluidRow(plotOutput('genes_plot'),
                   fluidRow(uiOutput("download_genes_plt_button")),
                   fluidRow(
                     column(3,offset=2,
                         sliderInput("genes",
                                     "Number of Genes Detected:",
                                     min=1,
                                     max=10000,
                                     value=1000)
                     ),
                     column(3,offset=2,
                         sliderInput("mitocontent",
                                     'Percent Mitochondria:',
                                     min=1,
                                     max=100,
                                     value=40)
                     )
                   )
                   )
  )),
  
  navbarMenu('Maps',
             tabPanel('PCA',
                      fluidRow(
                        column(12,
                               p(h3('principal component analysis'), "PCA projections of sample abundances onto any pair of components.")
                        ),
                        offset = 1),
                      fluidRow(
                        column(2,
                               selectInput('pc_x', label = 'x-axis PC: ', choices = 1:5,
                                           selected = 1)
                        ),
                        column(2,
                               selectInput('pc_y', label = 'y-axis PC: ', choices = 1:5,
                                           selected = 2)
                        )
                      ),
                      fluidRow(plotOutput('pca_plt'),click = "pca_click",
                               brush = brushOpts(
                                 id = "pca_brush")
                               ),
                      fluidRow(
                        div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                            downloadButton("download_pca_plt", "Download Plot"))),
                     
                       fluidRow(
                        column(12,
                               p(h3('loadings'), "observe contributions of transcripts to the principal component")
                        ),
                        offset = 1),
                      fluidRow(
                        column(3,
                               textInput('sample', label = 'transcript: ', value = ''
                               )
                        ),
                        column(2,
                               selectInput('pc_input', label = 'principal component: ', choices = 1:5,
                                           selected = 1)
                        ),
                        column(3,
                               selectInput('pc_count', label = 'number of PCs or transcripts: ', choices = 1:10,
                                           selected = 5))
                        ),
                      fluidRow(
                        column(12,
                               p(h3('variance explained'))
                        ),
                        offset = 1),
                      fluidRow(plotOutput('plt_pc_var')),
                      fluidRow(
                        column(12,
                               p(h3('loadings'))
                        ),
                        offset = 1),
                      fluidRow(plotOutput('plt_pc_loadings'))
             ),
             tabPanel('T-SNE',
                      fluidRow(
                        column(12,
                               p(h3('t-sne'),"t-Distributed Stochastic Neighbor Embedding")
                               ),
                        offset = 1),
                      fluidRow(
                        column(2,
                        selectInput('dim1', label = 'x-axis t-sne: ', choices = 1:3,
                                    selected = 1)
                        ),
                        column(2,
                             selectInput('dim2', label = 'y-axis t-sne: ', choices = 1:3,
                                         selected = 2)
                             
                      )
                      ),
                      fluidRow(plotOutput('tsne_plt'))
             )
  ),
  
  navbarMenu('DGE',
             tabPanel('Select-Groups',
                      fluidRow(
                        column(2,offset=1,
                               selectInput("cluster", "Choose a dataset:", 
                                           choices = c("PCA", "TSNE"))
                               ),
                        column(2,offset=1,
                               selectInput('dge_dim1',"x-axis:",choices = 1:5,selected=1)),
                        column(2,offset=1,
                               selectInput('dge_dim2','y-axis:',choices=1:5,selected =2))
                        ),
                      fluidRow(column(4,
                                      plotOutput('dge_plot1',brush = brushOpts(id="b1")
                                                 )
                                      ),
                               column(4,
                                      plotOutput('dge_plot2',brush=brushOpts(id='b2')))
                               ),
                      fluidRow(
                          h4('Cluster1',offset=1),
                          DT::dataTableOutput('cluster1')
                        ),
                      fluidRow(
                        h4('Cluster2',offset=1),
                        DT::dataTableOutput('cluster2')
                      ),
                      fluidRow(
                        h4('Top Differentially Expressed Genes',offset=1),
                        DT::dataTableOutput('dge')
                      )
                      ),
             tabPanel('MA-Plot',
                      fluidRow(plotOutput('maplot'))
                      ),
             tabPanel('Scatter-Plots'
                      ),
             tabPanel('Violin Plots'
                      )
             )

))

