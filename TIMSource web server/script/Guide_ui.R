Search_gene_ui = function(){
  tagList(
    br(),
    box(width = 12,
        h1("Parameter:"),
        div(style = "position:relative;left:35%;float:left;display:inline-block;width:30%;",
            div( 
              shinyInput_label_embed(tag = selectInput(inputId = "search_gene",label = "Please enter your gene of interest",choices = total_genes,selected = "TP53",width = "100%"),
                                     element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["gene_symbol_query",], placement = "top")),
              style = "display:inline-block;width:70%;"),
            div( actionButton(inputId = "search_single",label = "Submit",icon = icon("check-circle"),width = "100%"),style = "display:inline-block;width:25%;vertical-align: text-bottom;"),
        ),
        box(width = 12,
            strong('The Pancancer Investigation module allows users to easily analyze the relationship between gene mutations and immunotherapy as well as the tumor immune microenvironment (TIME) from a pan-cancer perspective. Simply enter the gene of interest in the selection box above and click the Submit button. The backend will perform a comprehensive analysis of the submitted gene across three major data sources, including meta-analysis and TIME-related analysis. To maintain the smoothness of the website, all analysis modules are turned off by default. Users can click the "+" on the right side of the corresponding module according to their needs, and the module will expand and initiate the analysis.')
            )
        # uiOutput(outputId = "gene_detial",style="display:inline-block;margin:50px")
        ),
    box(width = 12,
        tabsetPanel(id = "second_menus1",
          tabPanel(icon = icon('database'),
                   title = "ICB datasets",
                   br(),
                   column(width = 10,offset = 1,
                          div(style="border: 2px dashed grey;padding: 10px;margin-bottom:10px;padding-bottom: 0px;",
                            p("The Pancancer Investigation module provides a more clinically impactful approach from a pan-cancer perspective, allowing users to evaluate the associations of gene mutations and ICB outcomes (OS, PFS, response) through meta-analyses. In the ICB datasetss with transcriptional data, this module will also display the heatmap of the mean ssGSEA score difference of each signature between the mutation and wildtype tumors. The asterisks indicate statistical significance by Wilcoxon rank-sum test. ***p < 0.001, **p < 0.01, *p < 0.05. Users can click the '+' button to obtain the results according to their needs. To maintain the smoothness of the website, all analysis modules are turned off by default.")
                          )
                       ),
                   box(title = "Overall Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "OS_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Progress Free Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "PFS_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "RECIST",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RECIST_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Clinical Benefit",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RESPONSE_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "ICB datasets -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_single_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "ICB datasets -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_single_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
                   ),
          tabPanel(icon = icon('database'),
                   title = "TCGA datasets",
                   br(),
                   column(width = 10,offset = 1,
                          div(style="border: 2px dashed grey;padding: 10px;margin-bottom:10px;padding-bottom: 0px;",
                              p("The Pancancer Investigation module allows users to comprehensively evaluate the associations of gene mutations with immune features (immune infiltration and immune-related signatures) in each cancer type. It will display the heatmap of the mean ssGSEA score difference of each signature between the mutation and wildtype tumors. The asterisks indicate statistical significance by Wilcoxon rank-sum test. ***p < 0.001, **p < 0.01, *p < 0.05. Users can click the '+' button to obtain the results according to their needs.")
                          )
                   ),
                   box(title = "TCGA -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_single_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "TCGA -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_single_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
                   ),
          tabPanel(icon = icon('database'),
                   title = "CPTAC datasets",
                   br(),
                   column(width = 10,offset = 1,
                          div(style="border: 2px dashed grey;padding: 10px;margin-bottom:10px;padding-bottom: 0px;",
                              p("The Pancancer Investigation module allows users to comprehensively evaluate the associations of gene mutations with immune features (immune infiltration and immune-related signatures) in each cancer type. It will display the heatmap of the mean ssGSEA score difference of each signature between the mutation and wildtype tumors at both transcriptomic and proteomic levels. The asterisks indicate statistical significance by Wilcoxon rank-sum test. ***p < 0.001, **p < 0.01, *p < 0.05. Users can click the '+' button to obtain the results according to their needs.")
                          )
                   ),
                   box(title = "CPTAC -- Immune Infiltration(RNA)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_rna_single_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "CPTAC -- Immune Infiltration(Protein)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_protein_single_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "CPTAC -- Immune Signature(RNA)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_rna_single_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "CPTAC -- Immune Signature(Protein)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_protein_single_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
          )
          )
        ),

    div(style="height:200px;display: inline-block;width: 20px;")
  )
}

Search_pathway_ui = function(){
  tagList(
    br(),
    box(width = 12,
        h1("Parameter:"),
        div(
          div(style = "position:relative;left:15%;float:left;display:inline-block;width:33%;",
              div( 
                shinyInput_label_embed(tag = textInput(inputId = "search_pathway_kw",label = "Step 1: Input a keyword about your pathway of interest",width = "100%"),
                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["keyword",], placement = "top")),
                style = "display:inline-block;width:70%;"),
              div( actionButton(inputId = "search_pm_kw",label = "confirm",icon = icon("check-circle"),width = "100%"),style = "display:inline-block;width:25%")
          ),
          div(style = "position:relative;right:15%;float:right;display:inline-block;width:33%;",
              div( 
                shinyInput_label_embed(tag = pickerInput(inputId = "search_pathway",label = "Step2: Choose a Pathway",choices = NULL,multiple = F,width = "100%",options = list(title  = "Please input a Pathway")),
                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["pathway",], placement = "top")),
                style = "display:inline-block;width:70%;"),
              div( actionButton(inputId = "search_pm",label = "Submit",icon = icon("check-circle"),width = "100%"),style = "display:inline-block;width:25%")
            )
          ),
        box(width = 12,
            strong('The Pancancer Investigation module allows users to easily and quickly explore the associations between mutations in pathways or functions of interest and immunotherapy, as well as the tumor immune microenvironment (TIME) from a pan-cancer viewpoint. To begin, users should enter the keywords of the desired pathway or function into the left-side input box and submit. Subsequently, the system searches the MSigDB database for pathways or functions containing the keyword and presents them in a selection box on the right. Users can then select their pathway or function of interest from the right-side selection box and click the "Submit" button. The system carries out a comprehensive analysis using three primary data sources, including a meta-analysis and TIME-related analysis. To optimize the performance of website, all analysis modules are closed by default; however, users can expand their desired module by clicking the "+" symbol on the right side to initiate the analysis.')
        )
        ),
    
    box(width = 12,
        tabsetPanel(id = "second_menus1",
          tabPanel(icon = icon('database'),
                   title = "ICB datasets",
                   br(),
                   column(width = 10,offset = 1,
                          div(style="border: 2px dashed grey;padding: 10px;margin-bottom:10px;padding-bottom: 0px;",
                              p("The Pancancer Investigation module provides a more clinically impactful approach from a pan-cancer perspective, allowing users to evaluate the associations of altered pathway and ICB outcomes (OS, PFS, response) through meta-analyses. Users should first type a keyword of interests and select a pathway in the beside column. In the ICB datasetss with transcriptional data, this module will also display the heatmap of the mean ssGSEA score difference of each signature between the mutation and wildtype tumors. The asterisks indicate statistical significance by Wilcoxon rank-sum test. ***p < 0.001, **p < 0.01, *p < 0.05. Users can click the '+' button to obtain the results according to their needs.")
                          )
                   ),
                   box(title = "Overall Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "OS_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Progress Free Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "PFS_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "RECIST",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RECIST_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Clinical Benefit",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RESPONSE_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "ICB datasets -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_pm_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "ICB datasets -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_pm_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
          ),
          tabPanel(icon = icon('database'),
                   title = "TCGA datasets",
                   br(),
                   column(width = 10,offset = 1,
                          div(style="border: 2px dashed grey;padding: 10px;margin-bottom:10px;padding-bottom: 0px;",
                              p("The Pancancer Investigation module allows users to comprehensively evaluate the associations of altered pathway with immune features (immune infiltration and immune-related signatures) in each cancer type. It will display the heatmap of the mean ssGSEA score difference of each signature between the mutation and wildtype tumors. The asterisks indicate statistical significance by Wilcoxon rank-sum test. ***p < 0.001, **p < 0.01, *p < 0.05. Users can click the '+' button to obtain the results according to their needs.")
                          )
                   ),
                   box(title = "TCGA -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_pm_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "TCGA -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_pm_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
                   
          ),
          tabPanel(icon = icon('database'),
                   title = "CPTAC datasets",
                   br(),
                   column(width = 10,offset = 1,
                          div(style="border: 2px dashed grey;padding: 10px;margin-bottom:10px;padding-bottom: 0px;",
                              p("The Pancancer Investigation module allows users to comprehensively evaluate the associations of altered pathway with immune features (immune infiltration and immune-related signatures) in each cancer type. It will display the heatmap of the mean ssGSEA score difference of each signature between the mutation and wildtype tumors at both transcriptomic and proteomic levels. The asterisks indicate statistical significance by Wilcoxon rank-sum test. ***p < 0.001, **p < 0.01, *p < 0.05. Users can click the '+' button to obtain the results according to their needs.")
                          )
                   ),
                   box(title = "CPTAC -- Immune Infiltration(RNA)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_rna_pm_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "CPTAC -- Immune Infiltration(Protein)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_protein_pm_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "CPTAC -- Immune Signature(RNA)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_rna_pm_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "CPTAC -- Immune Signature(Protein)",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "CPTAC_protein_pm_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
          )
        )
    ),
    div(style="height:200px;display: flex;width: 20px;")
  )
}