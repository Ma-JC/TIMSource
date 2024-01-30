CPTAC_ui = function(){

  tagList(
    
    
    br(),
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        
        h2("Parameters:",actionBttn(inputId = "tutorial_cptac_single",label = "Guide",style = "unite",color = "danger",size = "sm"),style="color:white"),
        div(id="Cancer_type2_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Cancer_type2",label = "Cancer Type",choices = c("2019_Colon","2020_BI","2020_CCRCC","2020_HNSCC","2020_LSCC","2020_LUAD","2020_OVA","2020_UCEC","2021_GBM","2021_PDAC"),selected = "2020_CCRCC"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["cancer type cptac",], placement = "top"))
            
            ),
        div(id="gene_symbol_CPTAC_tutorial",
            shinyInput_label_embed(tag = selectizeInput(inputId = "gene_symbol_CPTAC",label = "Gene_Symbol",choices = NULL,multiple = F,width = "100%",options = list(placeholder = "Please input a gene")),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["gene_symbol_tcga",], placement = "top"))
            
            ),
        div(id="Mut_type_cptac_single_tutorial",
            shinyInput_label_embed(tag = pickerInput(inputId = "Mut_type_cptac_single",label = "Variant_Classification",choices = TCGA_mutation_type,selected = TCGA_mutation_type,multiple = T,options = list(`actions-box` = TRUE)),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Variant_Classification_tcga",], placement = "top"))
            
            ),
        div(id="Wild_type_cptac_single_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Wild_type_cptac_single",label = "Wildtype groups",choices = c("Wildtype","Others"),selected = "Wildtype",multiple = F,width = "100%"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Wildtype groups",], placement = "top"))
            
            ),
        actionButton(inputId = "sec3",label = "submit",icon = icon("arrow-alt-circle-up",id = "button_cptac_single"),style="color: #033c73;"),
        hr(),
        p("Single Gene -- CPTAC database",style="text-align:center;font-size=20px;font-weight:bolder;"),
        p("In the Browse module, the Single Gene function is designed to help users explore the relationship between the selected gene mutations and immunotherapy efficacy and immune microenvironment in a chosen cancer type and cohort. Based on your research design, you can choose a cancer type or cohort that meets your requirements and perform detailed analysis, gaining access to the analysis results in the form of graphs and tables. For specific information on each cancer type and cohort, please see the 'About' module. The CPTAC database provides multi-omics data from over 1,000 patients across 10 different cancer types or cohorts, which can be accessed from the Proteomic Data Commons (PDC). The biggest difference between the CPTAC database and the Ref_ICI datasets and TCGA database is that the CPTAC database provides proteomic data along with paired genomic and transcriptomic data. This allows us to explore the relationship between gene mutations and the immune microenvironment at the protein level and even compare the differences between transcriptomic and proteomic analysis results."),
        width = 3,id = "sidebar_id3"
        
      ),
      mainPanel = mainPanel(width = 8,id = "main_id3",
                            tabsetPanel(id = "cptac_single_analysis",
                              
                              tabPanel(title = "Overview",icon = icon("spinner"),
                                       dashboardPage(
                                         dashboardHeader(disable = T),
                                         dashboardSidebar(disable = T),
                                         dashboardBody(
                                           column(width = 12,
                                                  br(),
                                                  h2("CPTAC database | Overview",style="color:#033c73;"),
                                                  hr(style="background-color:#033c73;height:1px")),
                                           column(width = 12,
                                                  box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                                           ),
                                          
                                           box(title = "Plot1",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                               plotOutput(outputId = "maf1_CPTAC",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_cptac_single_uidown1")
                                                      )
                                               ),
                                           box(title = "Plot2",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                               plotOutput(outputId = "maf2_CPTAC",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_cptac_single_uidown2")
                                                      )
                                               ),
                                           box(title = "Plot3",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                               plotOutput(outputId = "maf3_CPTAC",width = "100%",height = "1000px"),
                                               column(width = 12,style = "margin-top:20px",
                                                      uiOutput(outputId = "over_cptac_single_uidown3")
                                               )
                                               ),
                                           a(id="button2",href="#",style="position: fixed ! important;right:20px;bottom:210px","back top")
                                         )
                                       )),
                              tabPanel(title = "Mutation",icon = icon("dna"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA database | Mutational Landscape",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                                       ),
                                       column(width = 12,bsAlert("warning_cptac_single")),
                                       box(title = "LollipopPlot",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                           withSpinner(plotOutput(outputId = "maf4_CPTAC",width = "100%",height = "1000px"),type = 1),
                                           column(width = 12,
                                                  uiOutput(outputId = "mut_cptac_single_uidown1")
                                                  )
                                           ),
                                       box(title = "LollipopPlot",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           withSpinner(plotOutput(outputId = "maf5_CPTAC",width = "100%",height = "1000px"),type = 1),
                                           column(width = 12,
                                                  uiOutput(outputId = "mut_cptac_single_uidown2")
                                                  ),
                                       ),
                                       box(title = "Mutation information",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           withSpinner(reactableOutput(outputId = "compare_mutation_CPTAC",width = "100%"),type = 1),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "mut_cptac_single_uitabledown")
                                           )
                                           
                                       )
                                       
                              ),
                              tabPanel(title = "Infiltrating immune cells",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Infiltrating immune cells",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                                       ),
                                       column(width = 12,bsAlert("warning2_cptac_single")),
                                       box(title = "Immune infiltration (RNA)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                           column(withSpinner(plotlyOutput(outputId = "immune_infiltration1_CPTAC_rna",height = "650px",width = "100%"),type = 1),width = 12),
                                           column(width = 12,
                                                  uiOutput(outputId = "infr_cptac_single_uidown")
                                                  )
                                       ),
                                       box(title = "Immune infiltration (Protein)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_infiltration1_CPTAC_protein",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "infp_cptac_single_uidown")
                                                )
                                       ),
                                       box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                         column(width = 12,
                                                br(),
                                                
                                                shinyInput_label_embed(tag = selectInput(inputId = "immune_cell_type_CPTAC",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                                                switchInput(inputId = "outlier2_CPTAC",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                                style="background-color:aliceblue;"),
                                         fluidRow(column(withSpinner(plotOutput(outputId = "immune_infiltration2_CPTAC_rna",height = "650px"),type = 1),width = 6),
                                                  column(withSpinner(plotOutput(outputId = "immune_infiltration2_CPTAC_protein",height = "650px"),type = 1),width = 6))
                                       )
                                       ),
                              
                              tabPanel(title = "Immune-related signatures",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Immune-related signatures",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                                       ),
                                       column(width = 12,bsAlert("warning3_cptac_single")),
                                       box(title = "Immune signatures (RNA)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_signature1_CPTAC_rna",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "sigr_cptac_single_uidown")
                                                )                                       ),
                                       box(title = "Immune signatures (Protein)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_signature1_CPTAC_protein",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "sigp_cptac_single_uidown")
                                                )
                                       ),
                                       box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                         column(width = 12,
                                                br(),
                                                
                                                shinyInput_label_embed(tag = selectInput(inputId = "immune_signature_CPTAC",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                                                                                '18-gene IFN signature',
                                                                                                                                                                                'Gene expression profile',
                                                                                                                                                                                'Cytolytic activity',
                                                                                                                                                                                '13 T-cell signature',
                                                                                                                                                                                'Effective T cell score',
                                                                                                                                                                                'Immune checkpoint expression',
                                                                                                                                                                                'TLS')
                                                                                         ,selected = "6-gene IFN signature",width = "100%"),
                                                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["sig",], placement = "top")),
                                                switchInput(inputId = "outlier3_CPTAC",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                                style="background-color:aliceblue;"),
                                         fluidRow(column(withSpinner(plotOutput(outputId = "immune_signature2_CPTAC_rna",height = "650px"),type = 1),width = 6),
                                                  column(withSpinner(plotOutput(outputId = "immune_signature2_CPTAC_protein",height = "650px"),type = 1),width = 6))
                                       ),

                                       
                                       
                              ),
                              tabPanel(title = "Differential expression analysis",icon = icon("table"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Differential expression analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                                       ),
                                       column(width = 12,bsAlert("warning4_cptac_single")),
                                       box(width = 12,title = "Volcano Plot (RNA)",solidHeader = T,collapsible = T,collapsed = F,
                                           fluidRow(width = 12,
                                                    br(),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct_CPTAC_rna",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "FC_CPTAC_rna",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue_CPTAC_rna",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    style="background-color:aliceblue;"
                                           ),
                                           
                                           column(withSpinner(plotlyOutput(outputId = "volcano_CPTAC_rna",height = "650px",width = "100%"),type = 1),width = 12),
                                           ),
                                       
                                       box(width = 12,title = "Differential expression gene (RNA)",solidHeader = T,collapsible = T,collapsed = T,
                                           column(width = 12,
                                                  br(),
                                                  withSpinner(reactableOutput(outputId = "DEG_tab_CPTAC_rna",width = "100%"),type = 1),
                                                  
                                           ),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "DEGr_cptac_single_uitabledown")
                                           )
                                       ),

                                       box(width = 12,title = "Volcano Plot (Protein)",solidHeader = T,collapsible = T,
                                           fluidRow(width = 12,
                                                    br(),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct_CPTAC_protein",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "FC_CPTAC_protein",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue_CPTAC_protein",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    style="background-color:aliceblue;"
                                           ),
                                           
                                           column(withSpinner(plotlyOutput(outputId = "volcano_CPTAC_protein",height = "650px",width = "100%"),type = 1),width = 12)
                                           ),
                                       box(width = 12,title = "Differential expression gene (Protein)",solidHeader = T,collapsible = T,collapsed = T,
                                         column(width = 12,
                                                br(),
                                                withSpinner(reactableOutput(outputId = "DEG_tab_CPTAC_protein",width = "100%"),type = 1),
                                         ),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "DEGp_cptac_single_uitabledown")
                                         )
                                       )
                                       
                              ),
                              tabPanel(
                                title = "GSEA",icon = icon("chart-area"),
                                column(width = 12,
                                       br(),
                                       h2("CPTAC database | GSEA",style="color:#033c73;"),
                                       hr(style="background-color:#033c73;height:1px")),
                                column(width = 12,
                                       box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                                ),
                                column(width = 12,bsAlert("warning5_cptac_single")),
                                box(
                                  width = 12,title = "GSEA(RNA)",solidHeader = T,collapsible = T,collapsed = T,
                                  fluidRow(width = 12,
                                           br(),
                                           column(
                                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_cptac_rna",
                                                                                            label = "Database",
                                                                                            choiceNames = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            choiceValues = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            justified = TRUE,selected = "HALLMARK",
                                                                                            checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                                                    element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["database",], placement = "top"))
                                             
                                                  ,width = 12),
                                           style="background-color:aliceblue;"
                                  ),
                                  column(width = 12,
                                         br(),
                                         div(id = "cptac_single_rna_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_cptac_rna",width = "100%",height = "720px"),type = 1))
                                         ),
                                         br(),
                                         br(),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "GSEAr_cptac_single_uitabledown")
                                         )
                                         
                                  ),
                                box(
                                  width = 12,title = "GSEA(Protein)",solidHeader = T,collapsible = T,collapsed = T,
                                  fluidRow(width = 12,
                                           br(),
                                           column(
                                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_cptac_protein",
                                                                                            label = "Database",
                                                                                            choiceNames = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            choiceValues = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            justified = TRUE,selected = "HALLMARK",
                                                                                            checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                                                    element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["database",], placement = "top"))
                                             
                                                  ,width = 12),
                                           style="background-color:aliceblue;"
                                  ),
                                  column(width = 12,
                                         br(),
                                         div(id = "cptac_single_protein_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_cptac_protein",width = "100%",height = "720px"),type = 1)),
                                  ),
                                  br(),
                                  br(),
                                  column(width = 12,style = "margin-top:20px",
                                         uiOutput(outputId = "GSEAp_cptac_single_uitabledown")
                                  )
                                  
                                )
                              ),
                              tabPanel(title = "Survival analysis",icon = icon("chart-line"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Survival analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                                       ),
                                       column(width = 12,bsAlert("warning6_cptac_single")),
                                       box(width = 12,title = "Survival analysis",solidHeader = T,collapsible = T,
                                           column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_survival",height = 650),type = 1)),
                                           br(),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "sur_cptac_single_uidown")
                                           )
                                           )

                                       
                              )
                              
                              
                            )
      )
      
    )
    
    
  )

}