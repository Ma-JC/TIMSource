TCGA_subtype_ui <- function(){
  tagList(
    
    br(),
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        
        
        h2("Parameters:",actionBttn(inputId = "tutorial_tcga_subtype",label = "Guide",style = "unite",color = "danger",size = "sm"),style="color:white"),
        div(id = "Cancer_type_subtype_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Cancer_type_subtype",label = "Cancer Type",choices = c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'),selected = "KIRC"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["cancer type tcga",], placement = "top"))
            
            ),
        div(id = "Subtype_tcga_tutorial",
            shinyInput_label_embed(tag = selectizeInput(inputId = "Subtype_tcga",label = "Mutation Subtype",choices = NULL),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Mutation Subtype",], placement = "top"))
            
            ),
        div(id = "Mut_type_tcga_subtype_tutorial",
            shinyInput_label_embed(tag = pickerInput(inputId = "Mut_type_tcga_subtype",label = "Variant_Classification",choices = TCGA_mutation_type,multiple = T,selected = TCGA_mutation_type,options = list(`actions-box` = TRUE)),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Variant_Classification_tcga",], placement = "top"))
            
            ),
        div(id = "Wild_type_tcga_subtype_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Wild_type_tcga_subtype",label = "Wildtype groups",choices = c("Wildtype","Others"),selected = "Wildtype",multiple = F,width = "100%"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Wildtype groups",], placement = "top"))
            
            ),
        div(id = "gene_symbol_TCGA_subtype_tutorial",
            shinyInput_label_embed(tag = selectizeInput(inputId = "gene_symbol_TCGA_subtype",label = "Gene_Symbol",choices = NULL,multiple = F,width = "100%",options = list(placeholder = "Please input one or more genes")),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Gene_Symbol",], placement = "top"))
            
            ),
        actionButton(inputId = "sec2_subtype",label = "submit",icon = icon("arrow-alt-circle-up",id = "button_tcga_subtype"),style="color: #033c73;"),
        hr(),
        p("Mutation Subtype -- TCGA database",style="text-align:center;font-size=20px;font-weight:bolder;"),
        p("Currently, gene-based transcriptome studies or web servers often don't take into account the patient's genetic mutation status. A large number of studies have shown that patients with different mutation backgrounds (especially driver mutations) have significant differences in gene regulation patterns and intercellular interactions, and even some genes have opposing prognostic values. Therefore, in the Browse module, the Mutation Subtype provides different gene mutation backgrounds for users to further compare and explore the relationship between gene expression and the immune microenvironment (or patient prognosis). The TCGA database originates from the Pan-Cancer Atlas (PanCanAtlas) initiative. Although it does not provide information on immune therapy efficacy, the dataset offers multi-omics data and clinical information for over 10,000 patients across 33 different cancer types. With this dataset, we can analyze the relationship between gene mutations and the immune microenvironment in a wider range of cancer types, thereby validating and expanding the results from the Ref_ICI datasets."),
        width = 3,id = "sidebar_id2_subtype"
        
      ),
      mainPanel = mainPanel(width = 8,id = "main_id2_subtype",
                            tabsetPanel(id = "tcga_subtype_analysis",
                              tabPanel(title = "Infiltrating immune cells",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA database | Infiltrating immune cells",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                                       ),
                                       column(width = 12,bsAlert("warning_tcga_subtype")),
                                       box(title = "Immune infiltration",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotOutput(outputId = "immune_infiltration1_subtype",height = "400px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "inf_tcga_subtype_uidown")
                                               )
                                       ),
                                       box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                         column(width = 12,
                                                br(),
                                                
                                                shinyInput_label_embed(tag = selectInput(inputId = "immune_cell_type_subtype",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                                                switchInput(inputId = "outlier2_subtype",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                                style="background-color:aliceblue;"),
                                         column(withSpinner(plotlyOutput(outputId = "immune_infiltration2_subtype",height = "650px",width = "100%"),type = 1),width = 12)
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
                                       column(width = 12,bsAlert("warning2_tcga_subtype")),
                                       box(title = "Immune signatures",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotOutput(outputId = "immune_signature1_subtype",height = "400px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "sig_tcga_subtype_uidown")
                                                ),
                                       ),
                                      box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                        column(width = 12,
                                               br(),
                                               
                                               shinyInput_label_embed(tag = selectInput(inputId = "immune_signature_subtype",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                                                                                 '18-gene IFN signature',
                                                                                                                                                                                 'Gene expression profile',
                                                                                                                                                                                 'Cytolytic activity',
                                                                                                                                                                                 '13 T-cell signature',
                                                                                                                                                                                 'Effective T cell score',
                                                                                                                                                                                 'Immune checkpoint expression',
                                                                                                                                                                                 'TLS')
                                                                                        ,selected = "6-gene IFN signature",width = "100%"),
                                                                      element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["sig",], placement = "top")),
                                               switchInput(inputId = "outlier3_subtype",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                               style="background-color:aliceblue;"),
                                        column(withSpinner(plotlyOutput(outputId = "immune_signature2_subtype",height = "650px",width = "100%"),type = 1),width = 12)
                                      )
                                       
                                       
                                       
                              ),
                              tabPanel(title = "Survival analysis",icon = icon("chart-line"),
                                       column(width = 12,
                                              br(),
                                              h2("TCGA database | Survival analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                                       ),
                                       column(width = 12,bsAlert("warning5_tcga_subtype")),
                                       box(width = 12,title = "Survival analysis",solidHeader = T,collapsible = T,collapsed = F,
                                           column(width = 12,withSpinner(plotOutput(outputId = "TCGA_survival_subtype",height = 1110),type = 1)),
                                           column(width = 12,
                                                  uiOutput(outputId = "sur_tcga_subtype_uidown")
                                           )
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
                                       column(width = 12,bsAlert("warning3_tcga_subtype")),
                                       box(width = 12,title = "Volcano Plot",solidHeader = T,collapsible = T,
                                           fluidRow(width = 12,
                                                    br(),
                                                    column(
                                                      shinyInput_label_embed(tag = selectInput(inputId = "MT_OR_WT_tcga",label = "Mutation or Wildtype",choices = c("All","Mutation","Wildtype"),selectize = FALSE,selected = "All"),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["WORM",], placement = "top"))
                                                      
                                                      ,width = 3),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct_subtype",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                                      
                                                      ,width = 3),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "FC_subtype",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      ,width = 3),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue_subtype",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                                      
                                                      ,width = 3),
                                                    style="background-color:aliceblue;"
                                           ),
                                           
                                           column(withSpinner(plotlyOutput(outputId = "volcano_subtype",height = "650px",width = "100%"),type = 1),width = 12),

                                       ),
                                       box(width = 12,title = "Differential expression gene",solidHeader = T,collapsible = T,collapsed = T,
                                         column(width = 12,
                                                br(),
                                                withSpinner(reactableOutput(outputId = "DEG_tab_subtype",width = "100%"),type = 1)
                                         ),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "DEG_tcga_subtype_uitabledown")
                                         )
                                       )
                              ),
                              tabPanel(title = "GSEA",icon = icon("chart-area"),
                                column(width = 12,
                                       br(),
                                       h2("TCGA database | GSEA",style="color:#033c73;"),
                                       hr(style="background-color:#033c73;height:1px")),
                                column(width = 12,
                                       box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                                ),
                                column(width = 12,bsAlert("warning4_tcga_subtype")),
                                box(
                                  width = 12,title = "GSEA",solidHeader = T,collapsible = T,collapsed = T,
                                  fluidRow(width = 12,
                                           br(),
                                           column(
                                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_tcga_subtype",
                                                                                            label = "Database",
                                                                                            choiceNames = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            choiceValues = c("GO_BP","GO_CC","GO_MF","KEGG","REACTOME","HALLMARK"),
                                                                                            justified = TRUE,selected = "HALLMARK",
                                                                                            checkIcon = list(yes = icon("ok",lib = "glyphicon"))),
                                                                    element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["database",], placement = "top")),
                                                  width = 12),
                                           # column(actionButton(inputId = "sub_pathway_gsea_tcga",label = "confirm",icon = icon("check-circle")),width = 2),
                                           style="background-color:aliceblue;"
                                  ),
                                  column(width = 12,
                                         br(),
                                         div(id="tcga_subtype_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_tcga_subtype",width = "100%",height = "720px"),type = 1))
                                         # downloadButton('TCGA_gsea_subtype_tabdown',label = 'Download Table')
                                  ),
                                  br(),
                                  br(),
                                  column(width = 12,style = "margin-top:20px",
                                         uiOutput(outputId = "GSEA_tcga_subtype_uitabledown")
                                  )
                                  
                                )
                              )

                              
                              
                            )
      )
      
    )
    
    
  )
}