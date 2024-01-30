Search_gene_ui = function(){
  tagList(
    br(),
    box(width = 12,
        h1("Parameter:"),
        div(style = "position:relative;left:35%;float:left;display:inline-block;width:30%;",
            div( 
              shinyInput_label_embed(tag = pickerInput(inputId = "search_gene",label = "Please enter your gene of interest",choices = total_genes,selected = "TP53",width = "100%"),
                                     element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["gene_symbol_query",], placement = "top")),
              style = "display:inline-block;width:70%;"),
            div( actionButton(inputId = "search_single",label = "confirm",icon = icon("check-circle"),width = "100%"),style = "display:inline-block;width:25%"),
        ),
        uiOutput(outputId = "gene_detial",style="display:inline-block;margin:50px")
        ),
    box(width = 12,
        tabsetPanel(id = "second_menus1",
          tabPanel(icon = icon('database'),
                   title = "Ref_ICI datasets",
                   br(),
                   box(title = "Overall Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "OS_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Progress Free Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "PFS_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "RECIST",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RECIST_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Clinical Benefit",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RESPONSE_single_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Ref_ICI -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_single_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Ref_ICI -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_single_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
                   ),
          tabPanel(icon = icon('database'),
                   title = "TCGA database",
                   br(),
                   box(title = "TCGA -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_single_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "TCGA -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_single_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
                   ),
          tabPanel(icon = icon('database'),
                   title = "CPTAC database",
                   br(),
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
          div(style = "position:relative;left:15%;float:left;display:inline-block;width:25%;",
              div( 
                shinyInput_label_embed(tag = textInput(inputId = "search_pathway_kw",label = "Step 1: Input a keyword about your pathway of interest",width = "100%"),
                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["keyword",], placement = "top")),
                style = "display:inline-block;width:70%;"),
              div( actionButton(inputId = "search_pm_kw",label = "confirm",icon = icon("check-circle"),width = "100%"),style = "display:inline-block;width:25%")
          ),
          div(style = "position:relative;right:15%;float:right;display:inline-block;width:25%;",
              div( 
                shinyInput_label_embed(tag = pickerInput(inputId = "search_pathway",label = "Step2: Choose a Pathway",choices = NULL,multiple = F,width = "100%",options = list(title  = "Please input a Pathway")),
                                       element = bs_embed_tooltip(tag = shiny_iconlink(name = "question-circle"),title = tooltip_text["pathway",], placement = "top")),
                style = "display:inline-block;width:70%;"),
              div( actionButton(inputId = "search_pm",label = "Submit",icon = icon("check-circle"),width = "100%"),style = "display:inline-block;width:25%")
            )
          )
        ),
    
    box(width = 12,
        tabsetPanel(id = "second_menus1",
          tabPanel(icon = icon('database'),
                   title = "Ref_ICI datasets",
                   br(),
                   box(title = "Overall Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "OS_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Progress Free Survival",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "PFS_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "RECIST",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RECIST_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Clinical Benefit",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "RESPONSE_pm_forest",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Ref_ICI -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_pm_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "Ref_ICI -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "ref_pm_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
          ),
          tabPanel(icon = icon('database'),
                   title = "TCGA database",
                   br(),
                   box(title = "TCGA -- Immune Infiltration",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_pm_Immune_Infiltration_heatmap",width = "100%",height = "1000px"),type = 1))),
                   box(title = "TCGA -- Immune Signature",solidHeader = T,collapsible = T,width = 12,collapsed = T,
                       column(width = 12,withSpinner(plotOutput(outputId = "TCGA_pm_Immune_pathway_heatmap",width = "100%",height = "1000px"),type = 1)))
                   
          ),
          tabPanel(icon = icon('database'),
                   title = "CPTAC database",
                   br(),
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