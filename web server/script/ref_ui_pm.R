ref_ui_pm = function(){tagList(
  br(), 
  sidebarLayout(
    sidebarPanel(
      h2("Parameters:",actionBttn(inputId = "tutorial_ref_pm",label = "Guide",style = "unite",color = "danger",size = "sm"),style="color:white"),
      div(id = "dataset_pm_tutorial",
          shinyInput_label_embed(tag = selectInput(inputId = "dataset_pm",label = "Step1: Choose a dataset",choices = c(
                                                                                                                        "1. Samstein et al, Bladder Cancer (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Breast Cancer (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Colorectal Cancer (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Esophagogastric Cancer (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Glioma (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Head and Neck Cancer (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Melanoma (Mixed ICB)",
                                                                                                                        "1. Samstein et al, NSCLC (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Renal Cell Carcinoma (Mixed ICB)",
                                                                                                                        "1. Samstein et al, Pan-cancer (Mixed ICB)",
                                                                                                                        "2. Van Allen et al, Melanoma (Anti-CTLA4)",
                                                                                                                        "3. Miao et al, Bladder Cancer (Mixed ICB)",
                                                                                                                        "3. Miao et al, Melanoma (Mixed ICB)",
                                                                                                                        "3. Miao et al, NSCLC (Mixed ICB)",
                                                                                                                        "3. Miao et al, Microsatellite-stable solid tumors (Mixed ICB)",
                                                                                                                        "4. Janjigian et al, Esophagogastric Cancer (Mixed ICB)",
                                                                                                                        "5. Rizvi et al, NSCLC (Mixed ICB)",
                                                                                                                        "6. Pender, Pan-cancer (Mixed ICB)",
                                                                                                                        "7. Hellmann, et al. NSCLC (Mixed ICB)",
                                                                                                                        "8. Hugo, et al. Melanoma (Anti-PD1/PDL1)",
                                                                                                                        "9. Jiao, et al, Gastrointestinal cancer (Mixed ICB)",
                                                                                                                        "10. Snyder, et al. Melanoma (Anti-CTLA4)",
                                                                                                                        "11. Mariathasan, et al. Urothelial cancer (Anti-PD1/PDL1)",
                                                                                                                        "12. Braun, et al. Clear cell renal cell carcinoma (Anti-PD1/PDL1)",
                                                                                                                        "13. Motzer, et al, renal cell carcinoma (Anti-PD1/PDL1+Axitinib)",
                                                                                                                        "14. Miao et al, Clear cell renal cell carcinoma (Anti-PD1/PDL1)",
                                                                                                                        "15. Zhao et al. Glioblastoma(Anti-PD1/PDL1)",
                                                                                                                        "16. Rizvi et al.(2015) NSCLC(Anti-PD1/PDL1)",
                                                                                                                        "17. Harding et al. Hepatocellular carcinoma(Mixed ICB)",
                                                                                                                        "18. Anagnostou et al. Melanoma(Mixed ICB)",
                                                                                                                        "19. Anagnostou et al. NSCLC(Mixed ICB)",
                                                                                                                        "20. Riaz et al. Melanoma(Mixed ICB)",
                                                                                                                        "21. Liu et al. Melanoma(Mixed ICB)",
                                                                                                                        "22. Bai et al. Gastric cancer(Mixed ICB)",
                                                                                                                        "23. Lu et al.  Neuroendocrine neoplasms(Anti-PD1/PDL1)",
                                                                                                                        "24. Wang et al. Gastrointestinal cancer(Mixed ICB)"
                                                                                                                      ),selected = "1. Samstein et al, Pan-cancer (Mixed ICB)"),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["dataset",], placement = "top"))
          
          ),
      div(id = "Keyword_tutorial",
          shinyInput_label_embed(tag = textInput(inputId = "Keyword",label = "Step2: Confirm a keyword"),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["keyword",], placement = "top"))
          ,style = "display:inline-block;width:70%"),
      div(id = "SUB_WORD_tutorial",actionButton(inputId = "SUB_WORD",label = "confirm",icon = icon("check-circle")),style = "display:inline-block;width:25%"),
      # div(id = "gene_symbol_pm_tutorial",selectizeInput(inputId = "gene_symbol_pm",label = "Step3: Choose a Pathway",choices = NULL,multiple = F,width = "100%",options = list(placeholder = "Please input a Pathway"))),
      div(id = "gene_symbol_pm_tutorial",
          shinyInput_label_embed(tag = pickerInput(inputId = "gene_symbol_pm",label = "Step3: Choose a Pathway",choices = NULL,multiple = F,width = "100%",options = list(title  = "Please input a Pathway")),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["pathway",], placement = "top"))
          
          ),
      div(id = "Mut_type_ref_pm_tutorial",
          shinyInput_label_embed(tag = pickerInput(inputId = "Mut_type_ref_pm",label = "Variant_Classification",choices = NULL,multiple = T,options = list(`actions-box` = TRUE)),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Variant_Classification",], placement = "top"))
          
          ),
      div(id = "Wild_type_ref_pm_tutorial",
          shinyInput_label_embed(tag = selectInput(inputId = "Wild_type_ref_pm",label = "Wildtype groups",choices = c("Wildtype","Others"),selected = "Wildtype",multiple = F,width = "100%"),
                                 element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Wildtype groups",], placement = "top"))
          
          ),
      actionButton(inputId = "sec_pm",label = "submit",icon = icon("arrow-alt-circle-up",id = "button_ref_pm"),style="color: #033c73;"),
      hr(),
      p("Pathway Mutation -- Ref_ICI datasets",style="text-align:center;font-size=20px;font-weight:bolder;"),
      p("In the Browse module, Pathway Mutation is designed to help users explore the relationship between interest pathway mutations and immunotherapy efficacy and immune microenvironment in selected cancer types and cohorts. The 'Pathway' contains gene sets from over 30,000 functions and pathways, sourced from the MSigDB database, and you can search for related pathways and functions through keywords. According to your research design, you can select a cancer type or cohort that meets your requirements for detailed analysis, and you can obtain analysis results in the form of graphs and tables. For details of each cancer type and cohort, please see the 'About' module. Ref_ICI datasets is a collection of 24 datasets related to immune checkpoint inhibitor treatments that provides valuable information on immunotherapy efficacy and some of the data also provides RNA-seq data, allowing us to analyze the relationship between gene mutations and immunotherapy efficacy and the immune microenvironment."),
      width = 3,id = "sidebar_id1_pm"),
    mainPanel(width = 8,id = "main_id1_pm",
              ############
              tabsetPanel(id = "ref_pm_analysis",
                tabPanel(
                  title = "Survival Outcomes",icon = icon("chart-line"),
                  column(width = 12,
                         br(),
                         h2("Ref_ICI datasets | Survival Outcomes",style="color:#033c73;"),
                         hr(style="background-color:#033c73;height:1px")),
                  column(width = 12,
                         box(width = 12,p("The Ref_ICI datasets provide immunotherapy-related survival analysis, including overall survival and progression-free survival. We use the Kaplan-Meier method to estimate the survival curve and the log-rank test to compare the difference in the survival curves between the mutation group and the wildtype group. Note that not all datasets contain both OS and PFS, and genes or pathway with less than 3 patients in either the mutation or wildtype group will not be analyzed."))
                  ),
                  column(width = 12,bsAlert("warning_pm")),
                  box(width = 12,title = "survival Curve",solidHeader = T,collapsible = T,
                    column(width = 12,withSpinner(plotOutput(outputId = "Survival_pm",height = 650),type = 1)),
                    column(width = 12,style = "margin-top:20px",
                           uiOutput(outputId = "sur_pm_uidown")
                           )
                  )

                ),
                tabPanel(title = "Drugs Response",icon = icon("capsules"),
                         column(width = 12,
                                br(),
                                h2("Ref_ICI datasets | Drugs Response",style="color:#033c73;"),
                                hr(style="background-color:#033c73;height:1px")),
                         column(width = 12,
                                box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                         ),
                         column(width = 12,bsAlert("warning2_pm")), 
                         box(width = 12,title = "Drugs Response Plot",solidHeader = T,collapsible = T,
                           column(width = 12,withSpinner(plotOutput(outputId = "response_pm",height = 650),type = 1)),
                           column(width = 12,style = "margin-top:20px",
                                  uiOutput(outputId = "res_pm_uidown")
                                  )
                         )

                ),
                tabPanel(title = "Mutation",icon = icon("dna"),
                         column(width = 12,
                                br(),
                                h2("Ref_ICI datasets | Mutation",style="color:#033c73;"),
                                hr(style="background-color:#033c73;height:1px")),
                         column(width = 12,
                                box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                         ),
                         column(width = 12,bsAlert("warning3_pm")),
                         box(width = 12,title = "Mutation Plot",solidHeader = T,collapsible = T,
                             column(switchInput(inputId = "outline_pm",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),width = 12),
                             column(width = 12,withSpinner(plotOutput(outputId = "TMB_pm",height = 650),type = 1)),
                             column(width = 12,style = "margin-top:20px",
                                    uiOutput(outputId = "mut_pm_uidown")
                                    )
                             ),
                         box(width = 12,title = "Patient Table",solidHeader = T,collapsible = T,collapsed = T,
                             column(width = 12,withSpinner(reactableOutput(outputId = "SPLOT_pm",width = "100%"),type = 1)),
                             column(width = 12,style = "margin-top:20px",
                                    uiOutput(outputId = "mut_pm_uitabledown")
                                    )
                             )

                         
                ),
                tabPanel(title = "Differential expression analysis",icon = icon("table"),
                         column(width = 12,
                                br(),
                                h2("Ref_ICI datasets | Differential expression analysis",style="color:#033c73;"),
                                hr(style="background-color:#033c73;height:1px")),
                         column(width = 12,
                                box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                         ),
                         column(width = 12,bsAlert("warning4_pm")),
                         box(width = 12,title = "Volcano Plot",solidHeader = T,collapsible = T,
                             fluidRow(width = 12,
                                      br(),
                                      column(
                                        shinyInput_label_embed(tag = numericInput(inputId = "min.pct_ref_pm",label = "Minimum percentage",value = 0.25,min = 0.05,max = 0.95,step = 0.05),
                                                               element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Minimum percentage",], placement = "top"))
                                        
                                        ,width = 4),
                                      column(
                                        shinyInput_label_embed(tag = numericInput(inputId = "FC_ref_pm",label = "|log2(Fold change)|",value = 1,min = 0.1,max = Inf,step = 0.1),
                                                               element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Fold change",], placement = "top"))
                                        
                                        ,width = 4),
                                      column(
                                        shinyInput_label_embed(tag = numericInput(inputId = "pvalue_ref_pm",label = "Adjust p value",value = 0.05,min = 0,max = 0.1,step = 0.005),
                                                               element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["Adjust p value",], placement = "top"))
                                        
                                        ,width = 4),
                                      style="background-color:aliceblue;"
                             ),
                             
                             column(withSpinner(plotlyOutput(outputId = "volcano_ref_pm",height = "650px",width = "100%"),type = 1),width = 12)
                         ),
                         box(width = 12,title = "Differential genes table",solidHeader = T,collapsible = T,collapsed = T,
                           column(width = 12,
                                  br(),
                                  withSpinner(reactableOutput(outputId = "DEG_tab_ref_pm",width = "100%"),type = 1)
                           ),
                           column(width = 12,style = "margin-top:20px",
                                  uiOutput(outputId = "DEG_pm_uitabledown")
                                  )
                         )
                         

                         
                ),
                tabPanel(
                  title = "GSEA",icon = icon("chart-area"),
                   column(width = 12,
                          br(),
                          h2("Ref_ICI datasets | GSEA",style="color:#033c73;"),
                          hr(style="background-color:#033c73;height:1px")),
                   column(width = 12,
                          box(width = 12,p("ssGSEA evaluated 28 types of immune cell infiltration in patients with different cancers, and observed the relationship between gene mutations and immune infiltration."))
                   ),
                  column(width = 12,bsAlert("warning5_pm")),
                   box(
                     width = 12,title = "GSEA",solidHeader = T,collapsible = T,collapsed = T,
                     fluidRow(width = 12,
                              br(),
                              column(
                                shinyInput_label_embed(tag = radioGroupButtons(inputId = "pathway_gsea_ref_pm",
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
                            div(id="ref_pm_GSEA_loader",style="height:100%",withSpinner(reactableOutput(outputId = "GSEA_tab_ref_pm",width = "100%",height = "720px"),type = 1)),
                            
                     ),
                     br(),
                     column(width = 12,style = "margin-top:20px",
                            uiOutput(outputId = "GSEA_pm_uitabledown")
                            )
                   )
                ),
                tabPanel(
                  title = "Immune Analysis",icon = icon("chart-bar"),
                  column(width = 12,
                         br(),
                         h2("Ref_ICI datasets | Immune Analysis",style="color:#033c73;"),
                         hr(style="background-color:#033c73;height:1px")),
                  column(width = 12,
                         box(width = 12,p("Patients with different cancers were scored on 8 immune-related signatures and observed the relationship between gene mutations and immunity."))
                  ),
                  column(width = 12,bsAlert("warning6_pm")),
                  box(width = 12,title = "Immune Infiltration",solidHeader = T,collapsible = T,
                    column(withSpinner(plotlyOutput(outputId = "immune_infiltration_datasets_all_pm",height = "650px",width = "100%"),type = 1),width = 12),
                    column(width = 12,style = "margin-top:20px",
                           uiOutput(outputId = "inf_pm_uidown")
                           )
                  ),
                  box(width = 12,title = "Immune Signature",solidHeader = T,collapsible = T,
                    column(withSpinner(plotlyOutput(outputId = "immune_signature_datasets_all_pm",height = "650px",width = "100%"),type = 1),width = 12),
                    column(width = 12,style = "margin-top:20px",
                           uiOutput(outputId = "sig_pm_uidown")
                           )
                  ),
                  box(width = 12,title = "Deital",solidHeader = T,collapsible = T,collapsed = T,
                    column(width = 6,
                           br(),
                           shinyInput_label_embed(tag = selectInput(inputId = "immune_infiltration_datasets_input_pm",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                           
                           shinyInput_label_embed(tag = selectInput(inputId = "immune_signature_datasets_input_pm",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
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
                           switchInput(inputId = "outlier_datasets_pm",label = "remove boxplot outlier",value = TRUE,width = "100%",onLabel = "YES",offLabel = "NO",labelWidth = "150px",handleWidth = "50px"),
                           style="background-color:aliceblue;"),
                    fluidRow(column(withSpinner(plotOutput(outputId = "immune_infiltration_datasets_one_pm",height = "650px"),type = 1),width = 6),
                             column(withSpinner(plotOutput(outputId = "immune_signature_datasets_one_pm",height = "650px"),type = 1),width = 6))
                  )

                  
                  
                  
                )
                
              )
              
    )
    
    
  ))}