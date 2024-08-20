ref_ui = function(){tagList(
  br(), 
  sidebarLayout(
    sidebarPanel(
      h2("Parameters:",actionBttn(inputId = "tutorial_ref_single",label = "Guide",style = "unite",color = "danger",size = "sm"),style="color:white"),
      div(id="dataset_tutorial",
          shinyInput_label_embed(tag = selectInput(inputId = "dataset",label = "Choose a dataset",choices = c(dataset_name3 ),selected = "Samstein et al, Pan-cancer (Mixed ICB)"),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["dataset",], placement = "top")
                                 )
          
          ),
      
      div(id = "gene_symbol_tutorial",
          shinyInput_label_embed(tag = selectizeInput(inputId = "gene_symbol",label = "Gene_Symbol",choices = NULL,multiple = F,width = "100%",options = list(placeholder = "Please input a gene")),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["gene_symbol",], placement = "top"))
          ),
      div(id = "Mut_type_ref_single_tutorial",
          shinyInput_label_embed(tag = pickerInput(inputId = "Mut_type_ref_single",label = "Variant_Classification",choices = NULL,multiple = T,options = list(`actions-box` = TRUE)),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Variant_Classification",], placement = "top"))
          ),
      # selectizeInput(inputId = "Mut_type_ref_single",label = "Variant_Classification",choices = NULL,multiple = T,width = "100%",options = list(placeholder = "Please input one or more mutation types")),
      div(id = "Wild_type_ref_single_tutorial",
          shinyInput_label_embed(tag = selectInput(inputId = "Wild_type_ref_single",label = "Wildtype groups",choices = c("Wildtype","Others"),selected = "Wildtype",multiple = F,width = "100%"),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Wildtype groups",], placement = "top"))

          ),
      div(id = "Therapy_type_ref_single_tutorial",
          shinyInput_label_embed(tag = pickerInput(inputId = "Therapy_type_ref_single",label = "Immune_checkpoint_blockade",choices = NULL,multiple = T,options = list(`actions-box` = TRUE)),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Immune_checkpoint_blockade",], placement = "top"))
      ),
      actionButton(inputId = "sec",label = "submit",icon = icon(name = "arrow-alt-circle-up",id = "button_ref_single"),style="color: #033c73;"),
      hr(),
      p("Single Gene",style="text-align:center;font-size=20px;font-weight:bolder;"),
      p("In the Inquire module, the Single Gene function is designed to help users explore the association of one single gene mutation with ICB outcomes, immune infiltration, immune-related signatures in a chosen cancer type and cohort. Users can choose one cancer type or cohort, and the mutation types to perform detailed analysis. The immune environment analysis is presented in TCGA and CPTAC section, and some of the ICB datasetss (if the transcriptional data was available in the published clinical study). For detailed information of the cohort, please see the 'About' module."),
      width = 3,id = "sidebar_id1"),
  mainPanel(width = 8,id = "main_id1",
            ############
            tabsetPanel(id = "ref_single_analysis",
              tabPanel(
                title = "Survival Outcomes",icon = icon("chart-line"),
                column(width = 12,
                       br(),
                       h2("ICB datasets | Survival Outcomes",style="color:#033c73;"),
                       hr(style="background-color:#033c73;height:1px")),
                column(width = 12,
                       box(width = 12,p("The ICB datasets provide immunotherapy-related survival analysis, including overall survival (OS) and progression-free survival (PFS). Please note that not all datasets contain both OS and PFS. The log-rank test is used to compare the survival differences between the mutation group and the wildtype group. Furthermore, if less than three patients in either mutation or wildtype group, the results will not be analyzed. This same criterion is applicable to other modules as well."))
                ),
                # column(width = 12,textOutput(outputId = "warning",inline = F)),
                column(width = 12,bsAlert("warning")),
                box(width = 12,title = "survival Curve",solidHeader = T,collapsible = T,
                    column(width = 12,withSpinner(plotOutput(outputId = "Survival",height = 650),type = 1)),
                    column(width = 12,style = "margin-top:20px",
                           uiOutput(outputId = "sur_single_uidown")
                           )
                    )

              ),
              tabPanel(title = "Drugs Response",icon = icon("capsules"),
                       column(width = 12,
                              br(),
                              h2("ICB datasets | Drugs Response",style="color:#033c73;"),
                              hr(style="background-color:#033c73;height:1px")),
                       column(width = 12,
                              box(width = 12,p("The ICB datasets efficacy data based on radiological response as per Response Evaluation Criteria in Solid Tumors (RECIST) provided by the original authors. Some cohorts employed a definition of 'clinical benefit', including patients with complete response (CR), partial response (PR), or stable disease (SD) lasting a durable time. We followed and utilized these definitions from original papers. The chi-square or Fisher's tests are used to compare the efficacy differences between the mutation group and the wildtype group. Please note that ICB efficacy information may not be available in all datasets."))
                       ),
                       # textOutput(outputId = "warning2",inline = F),  
                       column(width = 12,bsAlert("warning2")),
                       box(width = 12,title = "Drugs Response Plot",solidHeader = T,collapsible = T,
                           column(width = 12,withSpinner(plotOutput(outputId = "response",height = 650),type = 1)),
                           column(width = 12,style = "margin-top:20px",
                                  uiOutput(outputId = "res_single_uidown")
                                  )
                           )
              ),
              tabPanel(title = "Mutation",icon = icon("dna"),
                       column(width = 12,
                              br(),
                              h2("ICB datasets | Mutation",style="color:#033c73;"),
                              hr(style="background-color:#033c73;height:1px")),
                       column(width = 12,
                              box(width = 12,p("For the selected gene mutation, this module displays detailed mutation information including tumor mutation burden, mutation type, and detailed mutation information for each patient. Please note that the mutation profile is obtained from the post-processed data of individual study, and the mutation types may vary across different studies."))
                       ),
                       column(width = 12,bsAlert("warning3")),
                       box(width = 12,title = "Mutation Plot",solidHeader = T,collapsible = T,
                           column(switchInput(inputId = "outline",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),width = 12),
                           column(width = 12,withSpinner(plotOutput(outputId = "TMB",height = 650),type = 1)),
                           column(width = 12,style = "margin-top:20px",
                                  uiOutput(outputId = "mut_single_uidown")
                                  )
                           ),
                       box(width = 12,title = "Patient Table",solidHeader = T,collapsible = T,collapsed = T,
                         column(width = 12,
                                withSpinner(reactableOutput(outputId = "SPLOT",width = "100%"),type = 1)
                                ),
                         br(),
                         column(width = 12,style = "margin-top:20px",
                                uiOutput(outputId = "mut_single_uitabledown")
                                )
                       )

                       
              ),
              tabPanel(
                       title = "Differential expression analysis",
                       icon = icon("table"),
                       column(width = 12,
                              br(),
                              h2("ICB datasets | Differential expression analysis",style="color:#033c73;"),
                              hr(style="background-color:#033c73;height:1px")),
                       column(width = 12,
                              box(width = 12,p("Identifying the functional impact of somatic mutations on the gene expression pattern is critical for precision oncology and drug discovery. We provide the differential expression gene (DEG) analysis utilizing the limma package. This function allows users to set the parameters to perform customized analysis, aiming to identify the significant DEGs between the mutation group and the wild-type group. We present the specific information in the form of a volcano plot and table. Please note that this function can be performed in TCGA, CPTAC and some ICB datasets containing both genomics and transcriptomics data."))
                       ),
                       column(width = 12,bsAlert("warning4")),
                       box(width = 12,title = "Volcano Plot",solidHeader = T,collapsible = T,
                           fluidRow(width = 12,
                                    br(),
                                    column(
                                      shinyInput_label_embed(tag = numericInput(inputId = "min.pct_ref_single",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                      ,width = 4),
                                    column(
                                      shinyInput_label_embed(tag = numericInput(inputId = "FC_ref_single",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                      ,width = 4),
                                    column(
                                      shinyInput_label_embed(tag = numericInput(inputId = "pvalue_ref_single",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                             element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                      
                                      ,width = 4),
                                    style="background-color:aliceblue;"
                           ),
                           column(withSpinner(plotlyOutput(outputId = "volcano_ref_single",height = "650px",width = "100%"),type = 1),width = 12)
                       ),
                       box(width = 12,title = "Differential genes table",solidHeader = T,collapsible = T,collapsed = T,
                         column(width = 12,
                                br(),
                                withSpinner(reactableOutput(outputId = "DEG_tab_ref_single",width = "100%"),type = 1)
                                ),
                         br(),
                         column(width = 12,style = "margin-top:20px",
                                uiOutput(outputId = "DEG_single_uitabledown")
                                )
                       )
              ),
              tabPanel(
                title = "GSEA",icon = icon("chart-area"),
                column(width = 12,
                       br(),
                       h2("ICB datasets | GSEA",style="color:#033c73;"),
                       hr(style="background-color:#033c73;height:1px")),
                column(width = 12,
                       box(width = 12,p("In order to further understand the impact of gene mutations on function and pathways, we provide GSEA analysis supported by the clusterProfiler package. Functions/pathways are obtained from the MsigDB database including GO, KEGG, REACTOME, and HALLMARK. Results are presented in the table, including pathway name, set size, enrichment Score, NES, p value, p adjust value, q values, and rank, leading edge parameters and core enrichment genes. Clicking on the interested pathways in the presented table will lead to the displaying of GSEA plot, which can be downloaded through selected size, resolution, and multiple format types. Please be patient as GSEA analysis is time-consuming."))
                ),
                column(width = 12,bsAlert("warning5")),
                box(
                  width = 12,title = "GSEA",solidHeader = T,collapsible = T,collapsed = T,
                  fluidRow(width = 12,
                           br(),
                           column(
                             shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_ref_single",
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
                         div(id="ref_single_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_ref_single",width = "100%",height = "720px"),type = 1))
                        ),
                  br(),
                  column(width = 12,style = "margin-top:20px",
                         uiOutput(outputId = "GSEA_single_uitabledown")
                         )
                  
                  # column(withSpinner(plotOutput(outputId = "GSEA_plot_ref_single",height = "850px",width = "100%"),type = 1),width = 12),
                  # column(width = 12,
                  #        div(selectInput(inputId = "ref_gsea_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
                  #        div(numericInput(inputId = "ref_gsea_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
                  #        div(numericInput(inputId = "ref_gsea_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
                  #        div(numericInput(inputId = "ref_gsea_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
                  #        div(downloadButton(outputId = "ref_gsea_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;"))
                  # 
                )
              ),
              tabPanel(
                title = "Immune Analysis",icon = icon("chart-bar"),
                column(width = 12,
                       br(),
                       h2("ICB datasets | Immune Analysis",style="color:#033c73;"),
                       hr(style="background-color:#033c73;height:1px")),
                column(width = 12,
                       box(width = 12,p("To understand how the gene mutation orchestrate the tumor microenvironment (TME), we use the ssGSEA algorithm to evaluate the immune cell infiltration and immune-related features. The Wilcoxon rank-sum test is used to assess the statistical significance of the differences between the mutation group and the wild-type group."))
                ),
                column(width = 12,bsAlert("warning6")),
                box(width = 12,title = "Immune Infiltration",solidHeader = T,collapsible = T,
                  column(withSpinner(plotlyOutput(outputId = "immune_infiltration_datasets_all",height = "650px",width = "100%"),type = 1),width = 12),
                  column(width = 12,style = "margin-top:20px",
                         uiOutput(outputId = "inf_single_uidown"),
                         ),
                ),
                box(width = 12,title = "Immune Signature",solidHeader = T,collapsible = T,
                  column(withSpinner(plotlyOutput(outputId = "immune_signature_datasets_all",height = "650px",width = "100%"),type = 1),width = 12),
                  column(width = 12,style = "margin-top:20px",
                         uiOutput(outputId = "sig_single_uidown"),
                         ),
                  
                ),
                box(width = 12,title = "Deital",solidHeader = T,collapsible = T,collapsed = T,
                    column(width = 6,
                           br(),
                           shinyInput_label_embed(tag = selectInput(inputId = "immune_infiltration_datasets_input",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                           style="background-color:aliceblue;"),
                    column(width = 6,
                           br(),
                           
                           shinyInput_label_embed(tag = selectInput(inputId = "immune_signature_datasets_input",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                                                                    '18-gene IFN signature',
                                                                                                                                                                    'Gene expression profile',
                                                                                                                                                                    'Cytolytic activity',
                                                                                                                                                                    '13 T-cell signature',
                                                                                                                                                                    'Effective T cell score',
                                                                                                                                                                    'Immune checkpoint expression',
                                                                                                                                                                    'TLS')
                                                                    ,selected = "6-gene IFN signature",width = "100%"),
                                                  element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["sig",], placement = "top")),
                           style="background-color:aliceblue;"),
                    column(width = 12,
                           br(),
                           switchInput(inputId = "outlier_datasets",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                           style="background-color:aliceblue;"),
                    fluidRow(column(withSpinner(plotOutput(outputId = "immune_infiltration_datasets_one",height = "650px"),type = 1),width = 6),
                             column(withSpinner(plotOutput(outputId = "immune_signature_datasets_one",height = "650px"),type = 1),width = 6))
                    )

                
                
                
              )
              
            )
            
  )
  
  
))}