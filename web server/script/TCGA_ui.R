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
        p("Single Gene -- TCGA database",style="text-align:center;font-size=20px;font-weight:bolder;"),
        p("In the Browse module, the Single Gene function is designed to help users explore the relationship between the selected gene mutations and immunotherapy efficacy and immune microenvironment in a chosen cancer type and cohort. Based on your research design, you can choose a cancer type or cohort that meets your requirements and perform detailed analysis, gaining access to the analysis results in the form of graphs and tables. For specific information on each cancer type and cohort, please see the 'About' module. The TCGA database originates from the Pan-Cancer Atlas (PanCanAtlas) initiative. Although it does not provide information on immune therapy efficacy, the dataset offers multi-omics data and clinical information for over 10,000 patients across 33 different cancer types. With this dataset, we can analyze the relationship between gene mutations and the immune microenvironment in a wider range of cancer types, thereby validating and expanding the results from the Ref_ICI datasets."),
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
                                                  h2("TCGA database | Overview",style="color:#033c73;"),
                                                  hr(style="background-color:#033c73;height:1px")),
                                           column(width = 12,
                                                  box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                                           ),
                                           box(title = "Plot1",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                               plotOutput(outputId = "maf1",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_tcga_single_uidown1")
                                                      )
                                              ),
                                           box(title = "Plot2",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                               plotOutput(outputId = "maf2",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_tcga_single_uidown2")
                                                      )
                                              ),
                                           box(title = "Plot3",solidHeader = T,collapsible = T,collapsed = T,width = 12,
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
                                              h2("TCGA database | Mutational Landscape",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                                       ),
                                       column(width = 12,bsAlert("warning_tcga_single")),
                                       box(title = "LollipopPlot",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                           withSpinner(plotOutput(outputId = "maf4",width = "100%",height = "1000px"),type = 1),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "mut_tcga_single_uidown1")
                                                  )
                                           ),
                                       box(title = "LollipopPlot",solidHeader = T,collapsible = T,collapsed = T,width = 12,
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
                                              h2("TCGA database | Infiltrating immune cells",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
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
                                              h2("TCGA database | Immune-related signatures",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
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
                                              h2("TCGA database | Differential expression analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
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
                                       h2("TCGA database | GSEA",style="color:#033c73;"),
                                       hr(style="background-color:#033c73;height:1px")),
                                column(width = 12,
                                       box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
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
                                              h2("TCGA database | Survival analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
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
