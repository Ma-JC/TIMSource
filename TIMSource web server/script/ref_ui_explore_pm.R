ref_ui_explore_pm = function(){
  
  tagList(
    br(),
    box(width = 12,
        h2("Comprehensive Analysis of reference datasets",style="text-align:center;margin-top:10px;margin-bottom:10px;color:#033c73"),
        reactableOutput(outputId = "ref_datasets_detial_pm")
        ),
    box(width = 12,
        tabsetPanel(
          id = "explore_ref_pm",
          tabPanel(title = "Survival outcomes",
                   br(),
                   h3("Overall Survival",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   column(width = 12,bsAlert("warning_explore_pm_OS")),
                   column(width = 12,withSpinner(reactableOutput(outputId = "OS_total_pm"),type = 1)),
                   br(),
                   h3("Progress Free Survival",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   column(width = 12,bsAlert("warning_explore_pm_PFS")),
                   column(width = 12,withSpinner(reactableOutput(outputId = "PFS_total_pm"),type = 1)),
                   div(style="height:200px;display: inline-block;width: 20px;")
          
          ),
          tabPanel(title = "Drugs response",
                   br(),
                   h3("RECIST",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   column(width = 12,bsAlert("warning_explore_pm_recist")),
                   column(width = 12,withSpinner(reactableOutput(outputId = "RECIST_pm"),type = 1)),
                   br(),
                   h3("Clinical benefit",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   column(width = 12,bsAlert("warning_explore_pm_res")),
                   column(width = 12,withSpinner(reactableOutput(outputId = "RESPONSE_pm"),type = 1)),
                   div(style="height:200px;display: inline-block;width: 20px;")
         ),
          tabPanel(title = "Immune-associated analysis",
                   br(),
                   h3("Infiltrating immune cells",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   column(width = 12,bsAlert("warning_explore_pm_inf")),
                   column(width = 12,selectInput(inputId = "cell_type_pm",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                   withSpinner(reactableOutput(outputId = "immune_infiltrating_pm"),type = 1)),
                   br(),
                   h3("Immune-related signatures",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   column(width = 12,bsAlert("warning_explore_pm_sig")),
                   column(width = 12,selectInput(inputId = "signature_pm",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                               '18-gene IFN signature',
                                                                                                                               'Gene expression profile',
                                                                                                                               'Cytolytic activity',
                                                                                                                               '13 T-cell signature',
                                                                                                                               'Effective T cell score',
                                                                                                                               'Immune checkpoint expression',
                                                                                                                               'TLS')
                               ,selected = "6-gene IFN signature",width = "100%"),
                   withSpinner(reactableOutput(outputId = "immune_signature_pm"),type = 1)),
                   div(style="height:200px;display: inline-block;width: 20px;")
         ),
          div(style="height:200px")
        )
        )

  )
  
}