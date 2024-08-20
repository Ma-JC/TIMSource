TCGA_ui = function(){
  tagList(
    br(),
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        
        h2("Parameters:",actionBttn(inputId = "tutorial_tcga_single",label = "Guide",style = "unite",color = "danger",size = "sm"),style="color:white"),
        div(id = "Cancer_type_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Cancer_type",label = "Cancer Type",choices = c('ACC','BLCA','BRCA','CESC','CHOL','COAD',
                                                                                                               'DLBC','ESCA','GBM','HNSC','KICH','KIRC',
                                                                                                               'KIRP','LAML','LGG','LIHC','LUAD','LUSC',
                                                                                                               'MESO','OV','PAAD','PCPG','PRAD','READ',
                                                                                                               'SARC','SKCM','STAD','TGCT','THCA','THYM',
                                                                                                               'UCEC','UCS','UVM'),selected = "KIRC"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["cancer type tcga",], placement = "top"))
            
            ),
        div(id = "gene_symbol_TCGA_tutorial",
            shinyInput_label_embed(tag = selectizeInput(inputId = "gene_symbol_TCGA",label = "Gene_Symbol",choices = NULL,multiple = F,width = "100%",options = list(placeholder = "Please input a gene")),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["gene_symbol_tcga",], placement = "top"))
            
            ),
        div(id = "Mut_type_tcga_single_tutorial",
            shinyInput_label_embed(tag = pickerInput(inputId = "Mut_type_tcga_single",label = "Variant_Classification",choices = TCGA_mutation_type,selected = TCGA_mutation_type,multiple = T,options = list(`actions-box` = TRUE)),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Variant_Classification_tcga",], placement = "top"))
            
            ),
        div(id = "Wild_type_tcga_single_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Wild_type_tcga_single",label = "Wildtype groups",choices = c("Wildtype","Others"),selected = "Wildtype",multiple = F,width = "100%"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Wildtype groups",], placement = "top"))
            
            ),
        actionButton(inputId = "sec2",label = "submit",icon = icon("arrow-alt-circle-up",id = "button_tcga_single"),style="color: #033c73;"),
        hr(),
        p("Single Gene",style="text-align:center;font-size=20px;font-weight:bolder;"),
        p("In the Inquire module, the Single Gene function is designed to help users explore the association of one single gene mutation with ICB outcomes, immune infiltration, immune-related signatures in a chosen cancer type and cohort. Users can choose one cancer type or cohort, and the mutation types to perform detailed analysis. The immune environment analysis is presented in TCGA and CPTAC section, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study). For detailed information of the cohort, please see the 'About' module."),
        width = 3,id = "sidebar_id2"
        
      ),
      mainPanel = mainPanel(width = 8,id = "main_id2",
                            tabsetPanel(id = "tcga_single_analysis",
                              
                              tabPanel(title = "Overview",icon = icon("spinner"),
                                       dashboardPage(
                                         dashboardHeader(disable = T),
                                         dashboardSidebar(disable = T),
                                         dashboardBody(
                                           column(width = 12,
                                                  br(),
                                                  h2("TCGA datasets | Overview",style="color:#033c73;"),
                                                  hr(style="background-color:#033c73;height:1px")),
                                           column(width = 12,
                                                  box(width = 12,p("This page is designed to provide the mutation landscape of the selected cancer type, which is supported by the maftools package."))
                                           ),
                                           box(title = "PlotmafSummary",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                               plotOutput(outputId = "maf1",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_tcga_single_uidown1")
                                                      )
                                              ),
                                           box(title = "Waterfall Plot",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                               plotOutput(outputId = "maf2",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_tcga_single_uidown2")
                                                      )
                                              ),
                                           box(title = "PlotTiTv",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                               plotOutput(outputId = "maf3",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_tcga_single_uidown3")
                                                      )
                                               ),
                                           a(id="button",href="#",style="position: fixed ! important;right:20px;bottom:210px","back top")
                                         )
                                       )
                                       ),
                              tabPanel(title = "Mutation",icon = icon("dna"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA datasets | Mutational Landscape",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("For the selected gene, we use a lollipopPlot to display the gene's mutation rate, mutation type, and mutation localization. Additionally, we utilize the mafcompare function to analyze the differentially mutated genes between the mutation group and the wildtype, allowing us to identify co-occurring and mutually exclusive mutated genes with the selected gene."))
                                       ),
                                       column(width = 12,bsAlert("warning_tcga_single")),
                                       box(title = "LollipopPlot",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                           withSpinner(plotOutput(outputId = "maf4",width = "100%",height = "1000px"),type = 1),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "mut_tcga_single_uidown1")
                                                  )
                                           ),
                                       box(title = "Mafcompare",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           withSpinner(plotOutput(outputId = "maf5",width = "100%",height = "1000px"),type = 1),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "mut_tcga_single_uidown2")
                                                  )
                                           ),
                                       box(title = "Mutation information",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           withSpinner(reactableOutput(outputId = "compare_mutation",width = "100%"),type = 1),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "mut_tcga_single_uitabledown")
                                                  )
                                         
                                       )
                                       
                                       ),
                              tabPanel(title = "Infiltrating immune cells",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA datasets | Infiltrating immune cells",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("To understand how the gene mutation orchestrate the immune infiltration, we have collected 28 gene sets of immune cells and used the ssGSEA algorithm to evaluate the immune cell infiltration. The Wilcoxon rank-sum test is used to assess the statistical significance of the differences between the mutation group and the wild-type group."))
                                       ),
                                       column(width = 12,bsAlert("warning2_tcga_single")),
                                       box(title = "Immune infiltration",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_infiltration1",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "inf_tcga_single_uidown")
                                                )
                                       ),
                                       box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           column(width = 12,
                                                  br(),
                                                  
                                                  shinyInput_label_embed(tag = selectInput(inputId = "immune_cell_type",label = h4("Immune cell type"),choices = c('Activated B cell',
                                                                                                                                                                   'Activated CD4 T cell',
                                                                                                                                                                   'Activated CD8 T cell',
                                                                                                                                                                   'Central memory CD4 T cell',
                                                                                                                                                                   'Central memory CD8 T cell',
                                                                                                                                                                   'Effector memeory CD4 T cell',
                                                                                                                                                                   'Effector memeory CD8 T cell',
                                                                                                                                                                   'Gamma delta T cell',
                                                                                                                                                                   'Immature  B cell',
                                                                                                                                                                   'Memory B cell',
                                                                                                                                                                   'Regulatory T cell',
                                                                                                                                                                   'T follicular helper cell',
                                                                                                                                                                   'Type 1 T helper cell',
                                                                                                                                                                   'Type 17 T helper cell',
                                                                                                                                                                   'Type 2 T helper cell',
                                                                                                                                                                   'Activated dendritic cell',
                                                                                                                                                                   'CD56bright natural killer cell',
                                                                                                                                                                   'CD56dim natural killer cell',
                                                                                                                                                                   'Eosinophil',
                                                                                                                                                                   'Immature dendritic cell',
                                                                                                                                                                   'Macrophage','Mast cell',
                                                                                                                                                                   'MDSC','Monocyte',
                                                                                                                                                                   'Natural killer cell',
                                                                                                                                                                   'Natural killer T cell',
                                                                                                                                                                   'Neutrophil',
                                                                                                                                                                   'Plasmacytoid dendritic cell')
                                                                                                                                                                   ,selected = "Activated CD8 T cell",width = "100%"),
                                                                         element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Inf",], placement = "top")),
                                                  switchInput(inputId = "outlier2",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                                  style="background-color:aliceblue;"),
                                           column(withSpinner(plotOutput(outputId = "immune_infiltration2",height = "650px"),type = 1),width = 12)
                                           )
                                       ),
                              
                              tabPanel(title = "Immune-related signatures",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA datasets | Immune-related signatures",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("To understand how the gene mutation orchestrate the tumor microenvironment (TME), we have collected 8 immune-related signatures and used the ssGSEA algorithm to evaluate the immune-related features. The Wilcoxon rank-sum test is used to assess the statistical significance of the differences between the mutation group and the wild-type group."))
                                       ),
                                       column(width = 12,bsAlert("warning3_tcga_single")),
                                       box(title = "Immune signatures",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_signature1",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "sig_tcga_single_uidown")
                                                )
                                       ),
                                       box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                         column(width = 12,
                                                br(),
                                                
                                                shinyInput_label_embed(tag = selectInput(inputId = "immune_signature",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                                                                          '18-gene IFN signature',
                                                                                                                                                                          'Gene expression profile',
                                                                                                                                                                          'Cytolytic activity',
                                                                                                                                                                          '13 T-cell signature',
                                                                                                                                                                          'Effective T cell score',
                                                                                                                                                                          'Immune checkpoint expression',
                                                                                                                                                                          'TLS')
                                                                                                                                                                          ,selected = "6-gene IFN signature",width = "100%"),
                                                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["sig",], placement = "top")),
                                                switchInput(inputId = "outlier3",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                                style="background-color:aliceblue;"),
                                         column(withSpinner(plotOutput(outputId = "immune_signature2",height = "650px"),type = 1),width = 12)
                                       )

                                       
                                       
                              ),
                              tabPanel(title = "Differential expression analysis",icon = icon("table"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA datasets | Differential expression analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Identifying the functional impact of somatic mutations on the gene expression pattern is critical for precision oncology and drug discovery. We provide the differential expression gene (DEG) analysis utilizing the limma package. This function allows users to set the parameters to perform customized analysis, aiming to identify the significant DEGs between the mutation group and the wild-type group. We present the specific information in the form of a volcano plot and table. Please note that this function can be performed in TCGA, CPTAC and some ICB datasets containing both genomics and transcriptomics data."))
                                       ),
                                       column(width = 12,bsAlert("warning4_tcga_single")),
                                       box(width = 12,title = "Volcano Plot",solidHeader = T,collapsible = T,
                                           fluidRow(width = 12,
                                                    br(),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "FC",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      
                                                      ,width = 4),
                                                    style="background-color:aliceblue;"
                                           ),
                                           column(withSpinner(plotlyOutput(outputId = "volcano",height = "650px",width = "100%"),type = 1),width = 12)
                                           
                                       ),
                                       box(width = 12,title = "Differential expression gene",solidHeader = T,collapsible = T,collapsed = T,
                                         column(width = 12,
                                                br(),
                                                withSpinner(reactableOutput(outputId = "DEG_tab",width = "100%"),type = 1),
                                                
                                         ),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "DEG_tcga_single_uitabledown")
                                                )
                                       )
                                       
                              ),
                              tabPanel(
                                title = "GSEA",icon = icon("chart-area"),
                                column(width = 12,
                                       br(),
                                       h2("TCGA datasets | GSEA",style="color:#033c73;"),
                                       hr(style="background-color:#033c73;height:1px")),
                                column(width = 12,
                                       box(width = 12,p("In order to further understand the impact of gene mutations on function and pathways, we provide GSEA analysis supported by the clusterProfiler package. Functions/pathways are obtained from the MsigDB database including GO, KEGG, REACTOME, and HALLMARK. Results are presented in the table, including pathway name, set size, enrichment Score, NES, p value, p adjust value, q values, rank, and leading edge parameters. Clicking on the interested pathways in the presented table will lead to the displaying of GSEA plot, which can be downloaded through selected size, resolution, and multiple format types. Please be patient as GSEA analysis is time-consuming."))
                                ),
                                column(width = 12,bsAlert("warning5_tcga_single")),
                                box(
                                  width = 12,title = "GSEA",solidHeader = T,collapsible = T,collapsed = T,
                                  fluidRow(width = 12,
                                           br(),
                                           column(
                                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_tcga",
                                                                                            label = "Database",
                                                                                            choiceNames = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            choiceValues = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            justified = TRUE,selected = "HALLMARK",
                                                                                            checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                                                    element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["database",], placement = "top"))
                                             
                                                  ,width = 12),
                                           # column(actionButton(inputId = "sub_pathway_gsea_tcga",label = "confirm",icon = icon("check-circle")),width = 2),
                                           style="background-color:aliceblue;"
                                  ),
                                  column(width = 12,
                                         br(),
                                         div(id="tcga_single_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_tcga",width = "100%",height = "720px"),type = 1))
                                  ),
                                  br(),
                                  br(),
                                  column(width = 12,style = "margin-top:20px",
                                         uiOutput(outputId = "GSEA_tcga_single_uitabledown")
                                         )

                                ),
                                # column(withSpinner(plotOutput(outputId = "GSEA_plot_tcga",height = "850px",width = "100%"),type = 1),width = 12),
                                # column(width = 12,
                                #        div(selectInput(inputId = "TCGA_gsea_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
                                #        div(numericInput(inputId = "TCGA_gsea_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
                                #        div(numericInput(inputId = "TCGA_gsea_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
                                #        div(numericInput(inputId = "TCGA_gsea_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
                                #        div(downloadButton(outputId = "TCGA_gsea_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;"))
                              ),
                              tabPanel(title = "Survival analysis",icon = icon("chart-line"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA datasets | Survival analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("The TCGA PanCanAtlas provides overall survival (OS), disease-specific survival (DFS), and progression-free interval (PFI) data. The log-rank test is used to compare the survival differences between the mutation group and the wildtype group. Furthermore, if less than three patients in either mutation or wildtype group, the results will not be analyzed."))
                                       ),
                                       column(width = 12,bsAlert("warning6_tcga_single")),
                                       box(width = 12,title = "Survival analysis",solidHeader = T,collapsible = T,
                                         column(width = 12,withSpinner(plotOutput(outputId = "TCGA_survival",height = 650),type = 1)),
                                         br(),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "sur_tcga_single_uidown")
                                                )
                                       )
                                       
                              )
                              
                              
                            )
      )
      
    )
  )
}
