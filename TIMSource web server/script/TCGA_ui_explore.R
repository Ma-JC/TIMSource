TCGA_ui_explore = function(){
  
  tagList(
    br(),
    box(width = 12,
        h2("Comprehensive Analysis & Visualization of TCGA database",style="text-align:center;margin-top:10px;margin-bottom:10px;color:#033c73"),
        reactableOutput(outputId = "TCGA_detial_single")
    ),
    box(width = 12,
        tabsetPanel(
          id = "explore_TCGA_single",
          tabPanel(title = "Immune-associated analysis",
                   br(),
                   h3("Infiltrating immune cells",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   selectInput(inputId = "TCGA_cell_type_single",label = h4("Immune cell type"),choices = c('Activated B cell',
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
                   withSpinner(reactableOutput(outputId = "TCGA_immune_infiltrating_single"),type = 1),
                   br(),
                   h3("Immune-related signatures",style="color:#033c73;font-weight: bold;"),
                   hr(style="background-color:#033c73;height:1px"),
                   selectInput(inputId = "TCGA_signature_single",label = h4("Immune_related signatures"),choices = c('6-gene IFN signature',
                                                                                                                '18-gene IFN signature',
                                                                                                                'Gene expression profile',
                                                                                                                'Cytolytic activity',
                                                                                                                '13 T-cell signature',
                                                                                                                'Effective T cell score',
                                                                                                                'Immune checkpoint expression',
                                                                                                                'TLS')
                               ,selected = "6-gene IFN signature",width = "100%"),
                   withSpinner(reactableOutput(outputId = "TCGA_immune_signature_single"),type = 1)
          ),
          div(style="height:200px")
        )
    )
    
  )
  
}