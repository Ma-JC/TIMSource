CPTAC_ui_pm = function(){
  
  tagList(
    
    
    br(),
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        
        h2("Parameters:",actionBttn(inputId = "tutorial_cptac_pm",label = "Guide",style = "unite",color = "danger",size = "sm"),style="color:white"),
        div(id = "Cancer_type2_pm_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Cancer_type2_pm",label = "Step1: Choose a cancer Type",choices = c("2019_Colon","2020_BI","2020_CCRCC","2020_HNSCC","2020_LSCC","2020_LUAD","2020_OVA","2020_UCEC","2021_GBM","2021_PDAC"),selected = "2020_CCRCC"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["cancer type cptac",], placement = "top"))
            
            ),
        div(id = "Keyword_cptac_tutorial",
            shinyInput_label_embed(tag = textInput(inputId = "Keyword_cptac",label = "Step2: Confirm a keyword"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["keyword",], placement = "top"))
            
            ,style = "display:inline-block;width:70%"),
        div(id = "SUB_WORD_cptac_tutorial",actionButton(inputId = "SUB_WORD_cptac",label = "confirm",icon = icon("check-circle")),style = "display:inline-block;width:25%"),
        div(id = "gene_symbol_CPTAC_pm_tutorial",
            shinyInput_label_embed(tag = pickerInput(inputId = "gene_symbol_CPTAC_pm",label = "Step3: Choose a Pathway",choices = NULL,multiple = F,width = "100%",options = list(title = "Please input a pathway")),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["pathway",], placement = "top"))
            
            ),
        div(id = "Mut_type_cptac_pm_tutorial",
            shinyInput_label_embed(tag = pickerInput(inputId = "Mut_type_cptac_pm",label = "Variant_Classification",choices = TCGA_mutation_type,multiple = T,selected = TCGA_mutation_type,options = list(`actions-box` = TRUE)),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Variant_Classification_tcga",], placement = "top"))
            
            ),
        div(id = "Wild_type_cptac_pm_tutorial",
            shinyInput_label_embed(tag = selectInput(inputId = "Wild_type_cptac_pm",label = "Wildtype groups",choices = c("Wildtype","Others"),selected = "Wildtype",multiple = F,width = "100%"),
                                   element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Wildtype groups",], placement = "top"))
            
            ),
        actionButton(inputId = "sec3_pm",label = "submit",icon = icon("arrow-alt-circle-up",id = "button_cptac_pm"),style="color: #033c73;"),
        hr(),
        p("Pathway Mutation",style="text-align:center;font-size=20px;font-weight:bolder;"),
        p("In the Inquire module, 'Pathway Mutation' is designed to help users explore the association of altered pathway with ICB outcomes, immune infiltration, immune-related signatures in a chosen cancer type and cohort. The 'Pathway' contains gene sets from over 30,000 functions and pathways, sourced from the MSigDB database. By typing keywords, you can search for related pathways and select a pathway of interest to perform detailed analysis."),
        width = 3,
        id = "sidebar_id3_pm"
      ),
      mainPanel = mainPanel(width = 8,id = "main_id3_pm",
                            tabsetPanel(id = "cptac_pm_analysis",
                              
                              tabPanel(title = "Mutation",icon = icon("dna"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Mutational Landscape",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("For the selected pathway/process, we use a waterfall plot to display the mutation rate, mutation type, and mutation localization of each gene in the pathway. If there are more than 50 genes in the pathway, we only display the top 50 genes with the highest mutation frequency in the pathway. Additionally, we utilize the mafcompare function to analyze the differentially mutated genes between the mutation group and the wildtype, allowing us to identify co-occurring and mutually exclusive mutated genes with the selected gene."))
                                       ),
                                       column(width = 12,bsAlert("warning_cptac_pm")),
                                       box(title = "LollipopPlot",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                           withSpinner(plotOutput(outputId = "maf4_cptac_pm",width = "100%",height = "1500px"),type = 1),
                                           column(width = 12,
                                                  uiOutput(outputId = "mut_cptac_pm_uidown1")
                                                  )
                                           ),
                                       box(title = "Mafcompare",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           withSpinner(plotOutput(outputId = "maf5_cptac_pm",width = "100%",height = "1500px"),type = 1),
                                           column(width = 12,
                                                  uiOutput(outputId = "mut_cptac_pm_uidown2")
                                                  )

                                       ),
                                       box(title = "Mutation information",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                           withSpinner(reactableOutput(outputId = "compare_mutation_cptac_pm",width = "100%"),type = 1),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "mut_cptac_pm_uitabledown")
                                           )
                                           
                                       )
                                       
                              ),
                              tabPanel(title = "Infiltrating immune cells",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Infiltrating immune cells",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("To understand how the gene mutation orchestrate the immune infiltration, we have collected 28 gene sets of immune cells and used the ssGSEA algorithm to evaluate the immune cell infiltration. The Wilcoxon rank-sum test is used to assess the statistical significance of the differences between the mutation group and the wild-type group. Notably, CPTAC database contains the matched genomics, transcriptomics, and proteomics data, allowing immune infiltration analysis at both transcriptomic and proteomic level."))
                                       ),
                                       column(width = 12,bsAlert("warning2_cptac_pm")),
                                       box(title = "Immune infiltration (RNA)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_infiltration1_CPTAC_rna_pm",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "infr_cptac_pm_uidown")
                                                )
                                       ),
                                      box(title = "Immune infiltration (Protein)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                        column(withSpinner(plotlyOutput(outputId = "immune_infiltration1_CPTAC_protein_pm",height = "650px",width = "100%"),type = 1),width = 12),
                                        column(width = 12,
                                               uiOutput(outputId = "infp_cptac_pm_uidown")
                                               )
                                      ),
                                      box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                        column(width = 12,
                                               br(),
                                               
                                               shinyInput_label_embed(tag = selectInput(inputId = "immune_cell_type_CPTAC_pm",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                                               switchInput(inputId = "outlier2_CPTAC_pm",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                               style="background-color:aliceblue;"),
                                        fluidRow(column(withSpinner(plotOutput(outputId = "immune_infiltration2_CPTAC_rna_pm",height = "650px"),type = 1),width = 6),
                                                 column(withSpinner(plotOutput(outputId = "immune_infiltration2_CPTAC_protein_pm",height = "650px"),type = 1),width = 6))
                                      )
                                      
                                       
                                       
                              ),
                              
                              tabPanel(title = "Immune-related signatures",icon = icon("chart-bar"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Immune-related signatures",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("To understand how the gene mutation orchestrate the tumor microenvironment (TME), we have collected 8 immune-related signatures and used the ssGSEA algorithm to evaluate the immune-related features. The Wilcoxon rank-sum test is used to assess the statistical significance of the differences between the mutation group and the wild-type group. Notably, CPTAC database contains the matched genomics, transcriptomics, and proteomics data, allowing immune signature analysis at both transcriptomic and proteomic level."))
                                       ),
                                       column(width = 12,bsAlert("warning3_cptac_pm")),
                                       box(title = "Immune signatures (RNA)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_signature1_CPTAC_rna_pm",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "sigr_cptac_pm_uidown")
                                                )
                                       ),
                                       box(title = "Immune signatures (Protein)",solidHeader = T,collapsible = T,collapsed = F,width = 12,
                                         column(withSpinner(plotlyOutput(outputId = "immune_signature1_CPTAC_protein_pm",height = "650px",width = "100%"),type = 1),width = 12),
                                         column(width = 12,
                                                uiOutput(outputId = "sigp_cptac_pm_uidown")
                                               )
                                       ),
                                       box(title = "Detial",solidHeader = T,collapsible = T,collapsed = T,width = 12,
                                         column(width = 12,
                                                br(),
                                                
                                                shinyInput_label_embed(tag = selectInput(inputId = "immune_signature_CPTAC_pm",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                                                                                   '18-gene IFN signature',
                                                                                                                                                                                   'Gene expression profile',
                                                                                                                                                                                   'Cytolytic activity',
                                                                                                                                                                                   '13 T-cell signature',
                                                                                                                                                                                   'Effective T cell score',
                                                                                                                                                                                   'Immune checkpoint expression',
                                                                                                                                                                                   'TLS')
                                                                                                                                                                                   ,selected = "6-gene IFN signature",width = "100%"),
                                                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["sig",], placement = "top")),
                                                switchInput(inputId = "outlier3_CPTAC_pm",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                                                style="background-color:aliceblue;"),
                                         fluidRow(column(withSpinner(plotOutput(outputId = "immune_signature2_CPTAC_rna_pm",height = "650px"),type = 1),width = 6),
                                                  column(withSpinner(plotOutput(outputId = "immune_signature2_CPTAC_protein_pm",height = "650px"),type = 1),width = 6))
                                       )
                                       
                                      
                                       
                                       
                              ),
                              tabPanel(title = "Differential expression analysis",icon = icon("table"),
                                       column(width = 12,
                                              br(),
                                              h2("CPTAC database | Differential expression analysis",style="color:#033c73;"),
                                              hr(style="background-color:#033c73;height:1px")),
                                       column(width = 12,
                                              box(width = 12,p("Identifying the functional impact of somatic mutations on the gene expression pattern is critical for precision oncology and drug discovery. We provide the differential expression gene (DEG) and protein (DEP) analysis utilizing the limma package. This function allows users to set the parameters to perform customized analysis, aiming to identify the significant DEGs and DEPs between the mutation group and the wild-type group. We present the specific information in the form of a volcano plot and table."))
                                       ),
                                       column(width = 12,bsAlert("warning4_cptac_pm")),
                                       box(width = 12,title = "Volcano Plot (RNA)",solidHeader = T,collapsible = T,
                                           fluidRow(width = 12,
                                                    br(),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct_CPTAC_rna_pm",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "FC_CPTAC_rna_pm",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue_CPTAC_rna_pm",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    style="background-color:aliceblue;"
                                           ),
                                           
                                           column(withSpinner(plotlyOutput(outputId = "volcano_CPTAC_rna_pm",height = "650px",width = "100%"),type = 1),width = 12)
                                           ),
                                       box(width = 12,title = "Differential expression gene (RNA)",solidHeader = T,collapsible = T,collapsed = T,
                                         column(width = 12,
                                                br(),
                                                withSpinner(reactableOutput(outputId = "DEG_tab_CPTAC_rna_pm",width = "100%"),type = 1)
                                                
                                         ),
                                         column(width = 12,style = "margin-top:20px",
                                                uiOutput(outputId = "DEGr_cptac_pm_uitabledown")
                                         )

                                       ),

                                       box(width = 12,title = "Volcano Plot (Protein)",solidHeader = T,collapsible = T,collapsed = T,
                                           fluidRow(width = 12,
                                                    br(),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct_CPTAC_protein_pm",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "FC_CPTAC_protein_pm",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    column(
                                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue_CPTAC_protein_pm",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                                      
                                                      ,width = 4),
                                                    style="background-color:aliceblue;"
                                           ),
                                           
                                           column(withSpinner(plotlyOutput(outputId = "volcano_CPTAC_protein_pm",height = "650px",width = "100%"),type = 1),width = 12)

                                           ),
                                       box(width = 12,title = "Differential expression gene (Protein)",solidHeader = T,collapsible = T,collapsed = T,
                                           column(width = 12,
                                                  br(),
                                                  withSpinner(reactableOutput(outputId = "DEG_tab_CPTAC_protein_pm",width = "100%"),type = 1)
                                                  
                                           ),
                                           column(width = 12,style = "margin-top:20px",
                                                  uiOutput(outputId = "DEGp_cptac_pm_uitabledown")
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
                                       box(width = 12,p("In order to further understand the impact of gene mutations on function and pathways, we provide GSEA analysis supported by the clusterProfiler package. Functions/pathways are obtained from the MsigDB database including GO, KEGG, REACTOME, and HALLMARK. Results in transcriptomic and proteomic levels are presented in the table, including pathway name, set size, enrichment Score, NES, p value, p adjust value, q values, rank, and leading edge parameters. Clicking on the interested pathways in the presented table will lead to the displaying of GSEA plot, which can be downloaded through selected size, resolution, and multiple format types. Please be patient as GSEA analysis is time-consuming."))
                                ),
                                column(width = 12,bsAlert("warning5_cptac_pm")),
                                box(
                                  width = 12,title = "GSEA(RNA)",solidHeader = T,collapsible = T,collapsed = T,
                                  fluidRow(width = 12,
                                           br(),
                                           column(
                                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_cptac_rna_pm",
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
                                         div(id = "cptac_pm_rna_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_cptac_rna_pm",width = "100%",height = "720px"),type = 1))
                                         ),
                                  br(),
                                  br(),
                                  column(width = 12,style = "margin-top:20px",
                                         uiOutput(outputId = "GSEAr_cptac_pm_uitabledown")
                                  )
                                ),
                                box(
                                  width = 12,title = "GSEA(Protein)",solidHeader = T,collapsible = T,collapsed = T,
                                  fluidRow(width = 12,
                                           br(),
                                           column(
                                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_cptac_protein_pm",
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
                                         div(id = "cptac_pm_protein_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_cptac_protein_pm",width = "100%",height = "720px"),type = 1))
                                  ),
                                  br(),
                                  br(),
                                  column(width = 12,style = "margin-top:20px",
                                         uiOutput(outputId = "GSEAp_cptac_pm_uitabledown")
                                  )
                                  
                                )
                              )
                              # tabPanel(title = "Survival analysis",icon = icon("chart-line"),
                              #          column(width = 12,
                              #                 br(),
                              #                 h2("CPTAC database | Survival analysis",style="color:#033c73;"),
                              #                 hr(style="background-color:#033c73;height:1px")),
                              #          column(width = 12,
                              #                 box(width = 12,p("Due to the incomplete survival data obtained from CPTAC, the survival analysis may not fully align with the corresponding studies. For specific details, please refer to the original articles."))
                              #          ),
                              #          column(width = 12,bsAlert("warning6_cptac_pm")),
                              #          box(width = 12,title = "Survival analysis",solidHeader = T,collapsible = T,collapsed = F,
                              #            column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_survival_pm",height = 650),type = 1)),
                              #            column(width = 12,
                              #                   uiOutput(outputId = "sur_cptac_pm_uidown")
                              #                  )
                              #          )
                              #          )
                                       
                              
                              
                            )
      )
      
    )
    
    
  )
  
}