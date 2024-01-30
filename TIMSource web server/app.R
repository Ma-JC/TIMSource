# source("./script/loadData.R")
source("./script/loadFun.R")
source("./script/ref_ui.R")
source("./script/ref_ui_pm.R")
source("./script/ref_server.R")
source("./script/ref_server_pm.R")
source("./script/ref_ui_explore.R")
source("./script/ref_ui_explore_pm.R")
source("./script/TCGA_ui.R")
source("./script/TCGA_ui_pm.R")
source("./script/TCGA_server.R")
source("./script/TCGA_server_pm.R")
source("./script/TCGA_subtye_ui.R")
source("./script/TCGA_subtye_server.R")
source("./script/TCGA_ui_explore.R")
source("./script/TCGA_ui_explore_pm.R")
source("./script/CPTAC_ui.R")
source("./script/CPTAC_ui_pm.R")
source("./script/CPTAC_server.R")
source("./script/CPTAC_server_pm.R")
source("./script/CPTAC_subtype_ui.R")
source("./script/CPTAC_subtype_server.R")
source("./script/CPTAC_ui_explore.R")
source("./script/CPTAC_ui_explore_pm.R")
source("./script/Guide_ui.R")
source("./script/About_ui.R")
source("./script/All_tutorial.R")


ui <- fluidPage(id="fluid_id",title="SNVIO",theme = shinytheme("cerulean"),
                tags$head(tags$link(rel="shortcut icon",href="snvio.png")),
                use_cicerone(),
                shinyjs::useShinyjs(),
                use_bs_tooltip(),
                useWaiter(),
                navbarPage(title = p("Single Nucleotide Variant Immune Oncology (SNVIO)",style = "margin:0 0 0 25px"),
                           id = "navbar_id",
                           collapsible = T,
                           fluid = T,
                           inverse = T,
                            
                           #home#
                           tabPanel(title = icon("home"),
                                    value = "home",
                                    includeHTML(path = "www/shiny_home_v5.html")   
                           ),
                           
                           navbarMenu(title = "Broswe",
                                      icon = icon("arrows-alt"),
                                      tabPanel(title = "Single Gene",icon = icon("dna"),value = "Single gene",
                                               tabsetPanel(id = "second_menus1",
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "Ref_ICI datasets",
                                                             style = "padding:10px",
                                                             ref_ui()
                                                           ),
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "TCGA database",
                                                             style = "padding:10px",
                                                             TCGA_ui()
                                                           ),
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "CPTAC database",
                                                             style = "padding:10px",
                                                             CPTAC_ui()
                                                           )
                                                 )
                                               ),
                                      tabPanel(title = "Pathway Mutation",icon = icon("project-diagram"),
                                               tabsetPanel(id = "second_menus2",
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "Ref_ICI datasets",
                                                             style = "padding:10px",
                                                             ref_ui_pm()
                                                           ),
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "TCGA database",
                                                             style = "padding:10px",
                                                             TCGA_ui_pm()
                                                           ),
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "CPTAC database",
                                                             style = "padding:10px",
                                                             CPTAC_ui_pm()
                                                           )
                                                           
                                                           )
                                               ),
                                      tabPanel(title = "Mutation Subtype",icon = icon("network-wired"),
                                               tabsetPanel(id = "second_menus3",
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "TCGA database",
                                                             style = "padding:10px",
                                                             TCGA_subtype_ui()
                                                           ),
                                                           tabPanel(
                                                             icon = icon('database'),
                                                             title = "CPTAC database",
                                                             style = "padding:10px",
                                                             CPTAC_subtype_ui()
                                                           )
                                               )
                                               )
                                      
                                      
                          ),
                          #Explore#
                          navbarMenu(
                            title = "Explore",
                            icon = icon("table"),
                            tabPanel(title = "Single Gene",icon = icon("dna"),
                                     tabsetPanel(id = "second_menus1",
                                                 tabPanel(
                                                   icon = icon('database'),
                                                   title = "Ref_ICI datasets",
                                                   style = "padding:10px",
                                                   ref_ui_explore()


                                                 ),
                                                 tabPanel(
                                                   icon = icon('database'),
                                                   title = "TCGA database",
                                                   style = "padding:10px",
                                                   TCGA_ui_explore()
                                                 ),
                                                 tabPanel(
                                                   icon = icon('database'),
                                                   title = "CPTAC database",
                                                   style = "padding:10px",
                                                   CPTAC_ui_explore()

                                                 )
                                     )
                            ),
                            tabPanel(title = "Pathway Mutation",icon = icon("project-diagram"),
                                     tabsetPanel(id = "second_menus2",
                                                 tabPanel(
                                                   icon = icon('database'),
                                                   title = "Ref_ICI datasets",
                                                   style = "padding:10px",
                                                   ref_ui_explore_pm()

                                                 ),
                                                 tabPanel(
                                                   icon = icon('database'),
                                                   title = "TCGA database",
                                                   style = "padding:10px",
                                                   TCGA_ui_explore_pm()
                                                 ),
                                                 tabPanel(
                                                   icon = icon('database'),
                                                   title = "CPTAC database",
                                                   style = "padding:10px",
                                                   CPTAC_ui_explore_pm()

                                                 )

                                     )
                            )
                            # tabPanel(title = "Mutation Subtype",
                            #          tabsetPanel(id = "second_menus3",
                            #                      tabPanel(
                            #                        icon = icon('database'),
                            #                        title = "TCGA database",
                            #                        style = "padding:10px"
                            #
                            #                      ),
                            #                      tabPanel(
                            #                        icon = icon('database'),
                            #                        title = "CPTAC database",
                            #                        style = "padding:10px"
                            #
                            #                      )
                            #          )
                            # )



                          ),

						              #About#
						              navbarMenu(title = "Query",
						                         icon = icon("search"),
						                         tabPanel(title = "Single Gene",
						                                  icon = icon("dna"),
						                                  Search_gene_ui()
						                         ),
						                         tabPanel(title = "Pathway Mutation",
						                                  icon = icon("project-diagram"),
						                                  Search_pathway_ui()

						                         )
						              ),
						              tabPanel(title = "About",
						                       icon = icon("info-circle"),
						                        value = "About",
						                       style = "padding:10px",
						                       About_ui()
						              ),

						              div(style="height:200px"),
						              footer = footer(),
						              
						              style_snvio()
						              
                )
)

server <- function(input,output,session){
  


  ###################################HOME####################################################
  observeEvent(input$start,{
    updateTabsetPanel(session = session, inputId = "navbar_id", selected = "Single gene")
  })
  
  # observe({
  #   if(input$navbar_id == "Single gene"){
  #     message(isRunning())
  #   }
  # })

  ###################################Ref_datasets################################################
  ref = InputValidator$new()
  ref$add_rule("gene_symbol",sv_required())
  ref$add_rule("dataset",sv_required())
  ref$add_rule("Mut_type_ref_single",sv_required())
  ref$add_rule("Wild_type_ref_single",sv_required())
  ref$enable()
  
  ref2 = InputValidator$new()
  ref2$add_rule("min.pct_ref_single",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  ref2$add_rule("FC_ref_single",function(x){if(x <0){"warning: |log2FC| > 0"}})
  ref2$add_rule("pvalue_ref_single",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  ref2$enable()
  
  gene = eventReactive(input$sec,{req(ref$is_valid());input$gene_symbol})
  datasource  = eventReactive(input$sec,{req(ref$is_valid());dataset_name[input$dataset]})
  Mut_type_ref_single  = eventReactive(input$sec,{req(ref$is_valid());input$Mut_type_ref_single})
  Wild_type_ref_single  = eventReactive(input$sec,{req(ref$is_valid());input$Wild_type_ref_single})
  
  min.pct_ref_single  = eventReactive(input$min.pct_ref_single,{req(ref2$is_valid());input$min.pct_ref_single})
  FC_ref_single  = eventReactive(input$FC_ref_single,{req(ref2$is_valid());input$FC_ref_single})
  pvalue_ref_single  = eventReactive(input$pvalue_ref_single,{req(ref2$is_valid());input$pvalue_ref_single})
  
  observeEvent(input$dataset,{
    choices1 = gene_fre[[dataset_name[input$dataset]]]
    updateSelectizeInput(session,inputId = "gene_symbol",choices = choices1,selected = choices1[1],server = T)

    choices2 = unique(datasets_mu[[dataset_name[input$dataset]]]$Variant_Classification)
    updatePickerInput(session,inputId = "Mut_type_ref_single",choices = choices2,selected = choices2)

  })
  
  
  observeEvent(input$tutorial_ref_single,{
    guide$init()$start()
  })
  Ref_datasets_server(input,output,session,gene,datasource,datasets,datasets_mu,pathway_database,Mut_type_ref_single,Wild_type_ref_single,min.pct_ref_single,FC_ref_single,pvalue_ref_single)

  ###################################TCGA_datasets################################################
  tcga = InputValidator$new()
  tcga$add_rule("gene_symbol_TCGA",sv_required())
  tcga$add_rule("Cancer_type",sv_required())
  tcga$add_rule("Mut_type_tcga_single",sv_required())
  tcga$add_rule("Wild_type_tcga_single",sv_required())
  tcga$enable()
  
  tcga2 = InputValidator$new()
  tcga2$add_rule("min.pct",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  tcga2$add_rule("FC",function(x){if(x <1){if(x <0){"warning: |log2FC| > 0"}}})
  tcga2$add_rule("pvalue",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  tcga2$enable()
  
  
  gene_tcga = eventReactive(input$sec2,{req(tcga$is_valid());input$gene_symbol_TCGA})
  cancer_type  = eventReactive(input$sec2,{req(tcga$is_valid());input$Cancer_type})
  Mut_type_tcga_single  = eventReactive(input$sec2,{req(tcga$is_valid());input$Mut_type_tcga_single})
  Wild_type_tcga_single  = eventReactive(input$sec2,{req(tcga$is_valid());input$Wild_type_tcga_single})
  
  min.pct  = eventReactive(input$min.pct,{req(tcga2$is_valid());input$min.pct})
  FC  = eventReactive(input$FC,{req(tcga2$is_valid());input$FC})
  pvalue  = eventReactive(input$pvalue,{req(tcga2$is_valid());input$pvalue})
  
  observeEvent(input$Cancer_type,{
    #根据某一个输入来改变另一个输入的选项或数值，在shiny都是一类update开头的函数
    se = input$Cancer_type
    choices1 = TCGA[[se]]$maf@gene.summary$Hugo_Symbol
    updateSelectizeInput(session,inputId = "gene_symbol_TCGA",choices = choices1,selected = choices1[1],server = T)
  })

  observeEvent(input$tutorial_tcga_single,{
    guide_tcga_single$init()$start()
  })
  
  TCGA_server(input,output,session,gene_tcga,cancer_type,TCGA,pathway_database,Mut_type_tcga_single,Wild_type_tcga_single,min.pct,FC,pvalue)

  ###################################CPTAC_datasets################################################
  cptac = InputValidator$new()
  cptac$add_rule("gene_symbol_CPTAC",sv_required())
  cptac$add_rule("Cancer_type2",sv_required())
  cptac$add_rule("Mut_type_cptac_single",sv_required())
  cptac$add_rule("Wild_type_cptac_single",sv_required())
  cptac$enable()

  cptac2 = InputValidator$new()
  cptac2$add_rule("min.pct_CPTAC_rna",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  cptac2$add_rule("FC_CPTAC_rna",function(x){if(x <0){"warning: |log2FC| > 0"}})
  cptac2$add_rule("pvalue_CPTAC_rna",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  cptac2$enable()
  
  cptac3 = InputValidator$new()
  cptac3$add_rule("min.pct_CPTAC_protein",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  cptac3$add_rule("FC_CPTAC_protein",function(x){if(x <0){"warning: |log2FC| > 0"}})
  cptac3$add_rule("pvalue_CPTAC_protein",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  cptac3$enable()
  
  
  gene_cptac = eventReactive(input$sec3,{req(cptac$is_valid());input$gene_symbol_CPTAC})
  cancer_type2  = eventReactive(input$sec3,{req(cptac$is_valid());input$Cancer_type2})
  Mut_type_cptac_single  = eventReactive(input$sec3,{req(cptac$is_valid());input$Mut_type_cptac_single})
  Wild_type_cptac_single  = eventReactive(input$sec3,{req(cptac$is_valid());input$Wild_type_cptac_single})
  
  min.pct_CPTAC_rna  = eventReactive(input$min.pct_CPTAC_rna,{req(cptac2$is_valid());input$min.pct_CPTAC_rna})
  FC_CPTAC_rna  = eventReactive(input$FC_CPTAC_rna,{req(cptac2$is_valid());input$FC_CPTAC_rna})
  pvalue_CPTAC_rna  = eventReactive(input$pvalue_CPTAC_rna,{req(cptac2$is_valid());input$pvalue_CPTAC_rna})
  
  min.pct_CPTAC_protein  = eventReactive(input$min.pct_CPTAC_protein,{req(cptac3$is_valid());input$min.pct_CPTAC_protein})
  FC_CPTAC_protein  = eventReactive(input$FC_CPTAC_protein,{req(cptac3$is_valid());input$FC_CPTAC_protein})
  pvalue_CPTAC_protein  = eventReactive(input$pvalue_CPTAC_protein,{req(cptac3$is_valid());input$pvalue_CPTAC_protein})
  
  observeEvent(input$Cancer_type2,{
    #根据某一个输入来改变另一个输入的选项或数值，在shiny都是一类update开头的函数
    se = input$Cancer_type2
    choices1 = CPTAC[[se]]$maf@gene.summary$Hugo_Symbol
    updateSelectizeInput(session,inputId = "gene_symbol_CPTAC",choices = choices1,selected = choices1[1],server = T)

  })

  observeEvent(input$tutorial_cptac_single,{
    guide_cptac_single$init()$start()
  })
  
  CPTAC_server(input,output,session,gene_cptac,cancer_type2,CPTAC,pathway_database,Mut_type_cptac_single,Wild_type_cptac_single,
               min.pct_CPTAC_rna,FC_CPTAC_rna,pvalue_CPTAC_rna,min.pct_CPTAC_protein,FC_CPTAC_protein,pvalue_CPTAC_protein)


  ###################################Ref_datasets_pm################################################
  ref_pm = InputValidator$new()
  ref_pm$add_rule("gene_symbol_pm",sv_required())
  ref_pm$add_rule("dataset_pm",sv_required())
  ref_pm$add_rule("Mut_type_ref_pm",sv_required())
  ref_pm$add_rule("Wild_type_ref_pm",sv_required())
  # ref_pm$add_rule("Keyword",sv_required())
  ref_pm$enable()
  
  ref_pm2 = InputValidator$new()
  ref_pm2$add_rule("min.pct_ref_pm",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  ref_pm2$add_rule("FC_ref_pm",function(x){if(x <0){"warning: |log2FC| > 0"}})
  ref_pm2$add_rule("pvalue_ref_pm",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  ref_pm2$enable()
  
  

  gene_pm = eventReactive(input$sec_pm,{req(ref_pm$is_valid());input$gene_symbol_pm})
  datasource_pm  = eventReactive(input$sec_pm,{req(ref_pm$is_valid());dataset_name[input$dataset_pm]})
  Mut_type_ref_pm  = eventReactive(input$sec_pm,{req(ref_pm$is_valid());input$Mut_type_ref_pm})
  Wild_type_ref_pm  = eventReactive(input$sec_pm,{req(ref_pm$is_valid());input$Wild_type_ref_pm})
  
  min.pct_ref_pm  = eventReactive(input$min.pct_ref_pm,{req(ref_pm2$is_valid());input$min.pct_ref_pm})
  FC_ref_pm  = eventReactive(input$FC_ref_pm,{req(ref_pm2$is_valid());input$FC_ref_pm})
  pvalue_ref_pm  = eventReactive(input$pvalue_ref_pm,{req(ref_pm2$is_valid());input$pvalue_ref_pm})
  
  
  observeEvent(input$SUB_WORD,{
    target_pathway = names(pathway_list)[grepl(input$Keyword,names(pathway_list))]
    updatePickerInput(session,inputId = "gene_symbol_pm",choices = target_pathway,selected = target_pathway[1])
  })

  observeEvent(input$dataset_pm,{


      choices2 = unique(c(datasets_mu[[dataset_name[input$dataset_pm]]]$Variant_Classification))
      updatePickerInput(session,inputId = "Mut_type_ref_pm",choices = choices2,selected = choices2)


  })
  
  observeEvent(input$tutorial_ref_pm,{
    
    guide_pm$init()$start()
    
  })
  Ref_datasets_server_pm(input,output,session,gene_pm,datasource_pm,datasets,datasets_mu,Mut_type_ref_pm,Wild_type_ref_pm,min.pct_ref_pm,FC_ref_pm,pvalue_ref_pm)

  ###################################TCGA_datasets_pm################################################
  tcga_pm = InputValidator$new()
  tcga_pm$add_rule("gene_symbol_TCGA_pm",sv_required())
  tcga_pm$add_rule("Cancer_type_pm",sv_required())
  tcga_pm$add_rule("Mut_type_tcga_pm",sv_required())
  tcga_pm$add_rule("Wild_type_tcga_pm",sv_required())
  tcga_pm$enable()

  tcga_pm2 = InputValidator$new()
  tcga_pm2$add_rule("min.pct_pm",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  tcga_pm2$add_rule("FC_pm",function(x){if(x <0){"warning: |log2FC| > 0"}})
  tcga_pm2$add_rule("pvalue_pm",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  tcga_pm2$enable()
  
  gene_tcga_pm = eventReactive(input$sec2_pm,{req(tcga_pm$is_valid());input$gene_symbol_TCGA_pm})
  cancer_type_pm  = eventReactive(input$sec2_pm,{req(tcga_pm$is_valid());input$Cancer_type_pm})
  Mut_type_tcga_pm  = eventReactive(input$sec2_pm,{req(tcga_pm$is_valid());input$Mut_type_tcga_pm})
  Wild_type_tcga_pm  = eventReactive(input$sec2_pm,{req(tcga_pm$is_valid());input$Wild_type_tcga_pm})
  
  min.pct_pm  = eventReactive(input$min.pct_pm,{req(tcga_pm2$is_valid());input$min.pct_pm})
  FC_pm  = eventReactive(input$FC_pm,{req(tcga_pm2$is_valid());input$FC_pm})
  pvalue_pm  = eventReactive(input$pvalue_pm,{req(tcga_pm2$is_valid());input$pvalue_pm})
  
  observeEvent(input$SUB_WORD_tcga,{

    target_pathway = names(pathway_list)[grepl(input$Keyword_tcga,names(pathway_list))]
    updatePickerInput(session,inputId = "gene_symbol_TCGA_pm",choices = target_pathway,selected = target_pathway[1])
  })

  observeEvent(input$tutorial_tcga_pm,{
    
    guide_tcga_pm$init()$start()
    
  })
  
  TCGA_server_pm(input,output,session,gene_tcga_pm,cancer_type_pm,TCGA,pathway_database,Mut_type_tcga_pm,Wild_type_tcga_pm,min.pct_pm,FC_pm,pvalue_pm)
  ###################################CPTAC_datasets_pm################################################
  cptac_pm = InputValidator$new()
  cptac_pm$add_rule("gene_symbol_CPTAC_pm",sv_required())
  cptac_pm$add_rule("Cancer_type2_pm",sv_required())
  cptac_pm$add_rule("Mut_type_cptac_pm",sv_required())
  cptac_pm$add_rule("Wild_type_cptac_pm",sv_required())
  cptac_pm$enable()

  cptac_pm2 = InputValidator$new()
  cptac_pm2$add_rule("min.pct_CPTAC_rna_pm",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  cptac_pm2$add_rule("FC_CPTAC_rna_pm",function(x){if(x <0){"warning: |log2FC| > 0"}})
  cptac_pm2$add_rule("pvalue_CPTAC_rna_pm",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  cptac_pm2$enable()
  
  cptac_pm3 = InputValidator$new()
  cptac_pm3$add_rule("min.pct_CPTAC_protein_pm",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  cptac_pm3$add_rule("FC_CPTAC_protein_pm",function(x){if(x <0){"warning: |log2FC| > 0"}})
  cptac_pm3$add_rule("pvalue_CPTAC_protein_pm",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  cptac_pm3$enable()
  
  gene_cptac_pm = eventReactive(input$sec3_pm,{req(cptac_pm$is_valid());input$gene_symbol_CPTAC_pm})
  cancer_type2_pm  = eventReactive(input$sec3_pm,{req(cptac_pm$is_valid());input$Cancer_type2_pm})
  Mut_type_cptac_pm  = eventReactive(input$sec3_pm,{req(cptac_pm$is_valid());input$Mut_type_cptac_pm})
  Wild_type_cptac_pm  = eventReactive(input$sec3_pm,{req(cptac_pm$is_valid());input$Wild_type_cptac_pm})
  
  min.pct_CPTAC_rna_pm  = eventReactive(input$min.pct_CPTAC_rna_pm,{req(cptac_pm2$is_valid());input$min.pct_CPTAC_rna_pm})
  FC_CPTAC_rna_pm  = eventReactive(input$FC_CPTAC_rna_pm,{req(cptac_pm2$is_valid());input$FC_CPTAC_rna_pm})
  pvalue_CPTAC_rna_pm  = eventReactive(input$pvalue_CPTAC_rna_pm,{req(cptac_pm2$is_valid());input$pvalue_CPTAC_rna_pm})
  
  min.pct_CPTAC_protein_pm  = eventReactive(input$min.pct_CPTAC_protein_pm,{req(cptac_pm3$is_valid());input$min.pct_CPTAC_protein_pm})
  FC_CPTAC_protein_pm  = eventReactive(input$FC_CPTAC_protein_pm,{req(cptac_pm3$is_valid());input$FC_CPTAC_protein_pm})
  pvalue_CPTAC_protein_pm  = eventReactive(input$pvalue_CPTAC_protein_pm,{req(cptac_pm3$is_valid());input$pvalue_CPTAC_protein_pm})
  
  observeEvent(input$SUB_WORD_cptac,{

    target_pathway = names(pathway_list)[grepl(input$Keyword_cptac,names(pathway_list))]
    updatePickerInput(session,inputId = "gene_symbol_CPTAC_pm",choices = target_pathway,selected = target_pathway[1])
  })

  observeEvent(input$tutorial_cptac_pm,{
    
    guide_cptac_pm$init()$start()
    
  })
  
  CPTAC_server_pm(input,output,session,gene_cptac_pm,cancer_type2_pm,CPTAC,pathway_database,Mut_type_cptac_pm,Wild_type_cptac_pm,
                  min.pct_CPTAC_rna_pm,FC_CPTAC_rna_pm,pvalue_CPTAC_rna_pm,min.pct_CPTAC_protein_pm,FC_CPTAC_protein_pm,pvalue_CPTAC_protein_pm)


  ###################################TCGA_Subtype##################################################
  tcga_subtype = InputValidator$new()
  tcga_subtype$add_rule("Cancer_type_subtype",sv_required())
  tcga_subtype$add_rule("Subtype_tcga",sv_required())
  tcga_subtype$add_rule("gene_symbol_TCGA_subtype",sv_required())
  tcga_subtype$add_rule("Mut_type_tcga_subtype",sv_required())
  tcga_subtype$add_rule("Wild_type_tcga_subtype",sv_required())
  tcga_subtype$enable()

  tcga_subtype2 = InputValidator$new()
  tcga_subtype2$add_rule("MT_OR_WT_tcga",sv_required())
  tcga_subtype2$add_rule("min.pct_subtype",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  tcga_subtype2$add_rule("FC_subtype",function(x){if(x <0){"warning: |log2FC| > 0"}})
  tcga_subtype2$add_rule("pvalue_subtype",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  tcga_subtype2$enable()
  
  
  
  cancer_type_subtype  = eventReactive(input$sec2_subtype,{req(tcga_subtype$is_valid());input$Cancer_type_subtype})
  Mutation_subtype_tcga = eventReactive(input$sec2_subtype,{req(tcga_subtype$is_valid());input$Subtype_tcga})
  gene_tcga_subtype = eventReactive(input$sec2_subtype,{req(tcga_subtype$is_valid());input$gene_symbol_TCGA_subtype})
  Mut_type_tcga_subtype = eventReactive(input$sec2_subtype,{req(tcga_subtype$is_valid());input$Mut_type_tcga_subtype})
  Wild_type_tcga_subtype = eventReactive(input$sec2_subtype,{req(tcga_subtype$is_valid());input$Wild_type_tcga_subtype})
  
  MT_OR_WT_tcga = eventReactive(input$MT_OR_WT_tcga,{req(tcga_subtype2$is_valid());input$MT_OR_WT_tcga})
  min.pct_subtype = eventReactive(input$min.pct_subtype,{req(tcga_subtype2$is_valid());input$min.pct_subtype})
  FC_subtype = eventReactive(input$FC_subtype,{req(tcga_subtype2$is_valid());input$FC_subtype})
  pvalue_subtype = eventReactive(input$pvalue_subtype,{req(tcga_subtype2$is_valid());input$pvalue_subtype})
  
  
  observeEvent(input$Cancer_type_subtype,{
    #根据某一个输入来改变另一个输入的选项或数值，在shiny都是一类update开头的函数
    se = input$Cancer_type_subtype
    choices1 = TCGA[[se]]$maf@gene.summary$Hugo_Symbol
    choices2 = rownames(TCGA[[se]]$rna)
    updateSelectizeInput(session,inputId = "Subtype_tcga",choices = choices1,server = T)
    updateSelectizeInput(session,inputId = "gene_symbol_TCGA_subtype",choices = choices2,server = T)
  })

  observeEvent(input$tutorial_tcga_subtype,{
    
    guide_tcga_subtype$init()$start()
    
  })
  
  TCGA_subtype_server(input,output,session,cancer_type_subtype,Mutation_subtype_tcga,gene_tcga_subtype,TCGA,pathway_database,Mut_type_tcga_subtype,Wild_type_tcga_subtype,MT_OR_WT_tcga,min.pct_subtype,FC_subtype,pvalue_subtype)

  ###################################CPTAC_Subtype#################################################
  cptac_subtype = InputValidator$new()
  cptac_subtype$add_rule("Cancer_type2_subtype",sv_required())
  cptac_subtype$add_rule("Subtype_cptac",sv_required())
  cptac_subtype$add_rule("Mut_type_cptac_subtype",sv_required())
  cptac_subtype$add_rule("Wild_type_cptac_subtype",sv_required())
  cptac_subtype$add_rule("gene_symbol_CPTAC_subtype",sv_required())
  cptac_subtype$enable()
  
  cptac_subtype2 = InputValidator$new()
  cptac_subtype2$add_rule("MT_OR_WT_cptac_rna",sv_required())
  cptac_subtype2$add_rule("min.pct_CPTAC_rna_subtype",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  cptac_subtype2$add_rule("FC_CPTAC_rna_subtype",function(x){if(x <0){"warning: |log2FC| > 0"}})
  cptac_subtype2$add_rule("pvalue_CPTAC_rna_subtype",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  cptac_subtype2$enable()
  
  cptac_subtype3 = InputValidator$new()
  cptac_subtype3$add_rule("MT_OR_WT_cptac_protein",sv_required())
  cptac_subtype3$add_rule("min.pct_CPTAC_protein_subtype",function(x){if(x >0.95 | x< 0.05){"warning: 0.05 < min.pct < 0.95"}})
  cptac_subtype3$add_rule("FC_CPTAC_protein_subtype",function(x){if(x <0){"warning: |log2FC| > 0"}})
  cptac_subtype3$add_rule("pvalue_CPTAC_protein_subtype",function(x){if(x >0.1 | x< 0){"warning: 0 < pvalue < 0.1"}})
  cptac_subtype3$enable()

  cancer_type2_subtype  = eventReactive(input$sec3_subtype,{req(cptac_subtype$is_valid());input$Cancer_type2_subtype})
  Mutation_subtype_cptac = eventReactive(input$sec3_subtype,{req(cptac_subtype$is_valid());input$Subtype_cptac})
  gene_cptac_subtype = eventReactive(input$sec3_subtype,{req(cptac_subtype$is_valid());input$gene_symbol_CPTAC_subtype})
  Mut_type_cptac_subtype = eventReactive(input$sec3_subtype,{req(cptac_subtype$is_valid());input$Mut_type_cptac_subtype})
  Wild_type_cptac_subtype = eventReactive(input$sec3_subtype,{req(cptac_subtype$is_valid());input$Wild_type_cptac_subtype})
  
  MT_OR_WT_cptac_rna = eventReactive(input$MT_OR_WT_cptac_rna,{req(cptac_subtype2$is_valid());input$MT_OR_WT_cptac_rna})
  min.pct_CPTAC_rna_subtype = eventReactive(input$MT_OR_WT_cptac_rna,{req(cptac_subtype2$is_valid());input$min.pct_CPTAC_rna_subtype})
  FC_CPTAC_rna_subtype = eventReactive(input$FC_CPTAC_rna_subtype,{req(cptac_subtype2$is_valid());input$FC_CPTAC_rna_subtype})
  pvalue_CPTAC_rna_subtype = eventReactive(input$pvalue_CPTAC_rna_subtype,{req(cptac_subtype2$is_valid());input$pvalue_CPTAC_rna_subtype})
  
  MT_OR_WT_cptac_protein = eventReactive(input$MT_OR_WT_cptac_protein,{req(cptac_subtype3$is_valid());input$MT_OR_WT_cptac_protein})
  min.pct_CPTAC_protein_subtype = eventReactive(input$min.pct_CPTAC_protein_subtype,{req(cptac_subtype3$is_valid());input$min.pct_CPTAC_protein_subtype})
  FC_CPTAC_protein_subtype = eventReactive(input$FC_CPTAC_protein_subtype,{req(cptac_subtype3$is_valid());input$FC_CPTAC_protein_subtype})
  pvalue_CPTAC_protein_subtype = eventReactive(input$pvalue_CPTAC_protein_subtype,{req(cptac_subtype3$is_valid());input$pvalue_CPTAC_protein_subtype})
  
  observeEvent(input$Cancer_type2_subtype,{
    #根据某一个输入来改变另一个输入的选项或数值，在shiny都是一类update开头的函数
    se = input$Cancer_type2_subtype
    choices1 = CPTAC[[se]]$maf@gene.summary$Hugo_Symbol
    choices2 = intersect(rownames(CPTAC[[se]]$rna),rownames(CPTAC[[se]]$protein))
    updateSelectizeInput(session,inputId = "Subtype_cptac",choices = choices1,server = T)
    updateSelectizeInput(session,inputId = "gene_symbol_CPTAC_subtype",choices = choices2,server = T)
  })

  observeEvent(input$tutorial_cptac_subtype,{
    
    guide_cptac_subtype$init()$start()
    
  })

  CPTAC_subtype_server(input,output,session,cancer_type2_subtype,Mutation_subtype_cptac,gene_cptac_subtype,CPTAC,pathway_database,Mut_type_cptac_subtype,Wild_type_cptac_subtype,
                       MT_OR_WT_cptac_rna,min.pct_CPTAC_rna_subtype,FC_CPTAC_rna_subtype,pvalue_CPTAC_rna_subtype,MT_OR_WT_cptac_protein,min.pct_CPTAC_protein_subtype,FC_CPTAC_protein_subtype,pvalue_CPTAC_protein_subtype)



  ###################################Ref_datasets_Explore_single########################################
  output$ref_datasets_detial_single = renderReactable({
    datasets_overview$details = NA
    reactable(
      datasets_overview,
      rownames = F,
      searchable = TRUE,
      paginationType = "jump",
      resizable = TRUE,
      showPageSizeOptions = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      outlined = TRUE,
      striped = TRUE,
      defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                             sortNALast = TRUE,
                             align = "center"
                             # cell = function(value) format(value,digits = 2)
      ),
      columns = list(
        # Render a "show details" button in the last column of the table.
        # This button won't do anything by itself, but will trigger the custom
        # click action on the column.
        details = colDef(
          name = "",
          sortable = FALSE,
          cell = function() htmltools::tags$button("Show details")
        )
      ),
      onClick = JS("function(rowInfo, column) {
      if (column.id !== 'details') {
        return
      }
      if (window.Shiny) {
        Shiny.setInputValue('show_details_single', { index: rowInfo.index + 1 }, { priority: 'event' })
        }
                     }"
      )
    )
  })

  observeEvent(input$show_details_single$index,{

    ds = ref_name[input$show_details_single$index]

    if(ds %in% c("dataset5","dataset7","dataset13","dataset16")){
      closeAlert(session,"warning_explore_single_OS2")
      createAlert(session, "warning_explore_single_OS", "warning_explore_single_OS2", title = "Tips",
                  content = paste(ds,"didn't provide overall survival"), append = FALSE)
      closeAlert(session,"warning_explore_single_PFS2")
    }else if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                       "dataset1","dataset8","dataset10","dataset11")){
      closeAlert(session,"warning_explore_single_PFS2")
      createAlert(session, "warning_explore_single_PFS", "warning_explore_single_PFS2", title = "Tips",
                  content = paste(ds,"didn't provide progression-free survival"), append = FALSE)
      closeAlert(session,"warning_explore_single_OS2")
    }else{
      closeAlert(session,"warning_explore_single_OS2")
      closeAlert(session,"warning_explore_single_PFS2")
    }

  })
  output$OS_total_single = renderReactable({
    validate(need(!ref_name[input$show_details_single$index] %in% c("dataset5","dataset7","dataset13","dataset16"),message = FALSE))
    tmp = ref_total_OS_single[[ref_name[input$show_details_single$index]]]
    # tmp$Plot = NA
    reactable( tmp,
               searchable = TRUE,
               paginationType = "jump",
               resizable = TRUE,
               showPageSizeOptions = TRUE,
               #             onClick = "expand",
               highlight = TRUE,
               bordered = TRUE,
               outlined = TRUE,
               striped = TRUE,
               defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                      sortNALast = TRUE,
                                      align = "center"
                                      # cell = function(value) format(value,digits = 2)
               ),
               # columns = list(
               #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
               #   `Log_rank_test(OS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `Wald_test(OS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `HR(OS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
               #     list(color = color, fontWeight = "bold")}
               #   )
               # ),
               # onClick = JS(
               #   "function(rowInfo, column){
               #                      if (column.id !== 'Plot') {return}
               #                      window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
               #                      if (window.Shiny) {
               #                                          Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
               #                                        }
               #           }"
               # ),
               defaultSortOrder = "asc",
               defaultSorted = list(`Log_rank_test(OS)` = "asc", `HR(OS)` = "asc")

    )
  })
  output$PFS_total_single = renderReactable({
    validate(need(!ref_name[input$show_details_single$index] %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                                             "dataset1","dataset8","dataset10","dataset11"),message = FALSE))

    tmp = ref_total_PFS_single[[ref_name[input$show_details_single$index]]]
    # tmp$Plot = NA
    reactable( tmp,
               searchable = TRUE,
               paginationType = "jump",
               resizable = TRUE,
               showPageSizeOptions = TRUE,
               #             onClick = "expand",
               highlight = TRUE,
               bordered = TRUE,
               outlined = TRUE,
               striped = TRUE,
               defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                      sortNALast = TRUE,
                                      align = "center"
                                      # cell = function(value) format(value,digits = 2)
               ),
               # columns = list(
               #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
               #   `Log_rank_test(PFS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `Wald_test(PFS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `HR(PFS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
               #     list(color = color, fontWeight = "bold")}
               #   )
               # ),
               # onClick = JS(
               #   "function(rowInfo, column){
               #                      if (column.id !== 'Plot') {return}
               #                      window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
               #                      if (window.Shiny) {
               #                                          Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
               #                                        }
               #           }"
               # ),
               defaultSortOrder = "asc",
               defaultSorted = list(`Log_rank_test(PFS)` = "asc", `HR(PFS)` = "asc")

    )
  })


  observeEvent(input$show_details_single$index,{

    ds = ref_name[input$show_details_single$index]

    if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                 "dataset1","dataset13")){
      closeAlert(session,"warning_explore_single_recist2")
      closeAlert(session,"warning_explore_single_res2")
      createAlert(session, "warning_explore_single_recist", "warning_explore_single_recist2", title = "Warning",style = "danger",
                  content = paste(ds,"didn't provide data about RECIST"), append = FALSE)
      createAlert(session, "warning_explore_single_res", "warning_explore_single_res2", title = "Warning",style = "danger",
                  content = paste(ds,"didn't provide data about DCB"), append = FALSE)

    }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
      closeAlert(session,"warning_explore_single_recist2")
      createAlert(session, "warning_explore_single_recist", "warning_explore_single_recist2", title = "Tips",
                  content = paste(ds,"didn't provide data about RECIST"), append = FALSE)
      closeAlert(session,"warning_explore_single_res2")
    }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
      closeAlert(session,"warning_explore_single_res2")
      createAlert(session, "warning_explore_single_res", "warning_explore_single_res2", title = "Tips",
                  content = paste(ds,"didn't provide data about DCB"), append = FALSE)
      closeAlert(session,"warning_explore_single_recist2")
    }else{
      closeAlert(session,"warning_explore_single_recist2")
      closeAlert(session,"warning_explore_single_res2")
    }



  })
  output$RECIST_single = renderReactable({
    validate(need(!ref_name[input$show_details_single$index] %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                                             "dataset1","dataset13","dataset5","dataset6","dataset9","dataset10","dataset15","dataset19"),message = FALSE))
    tmp = ref_total_RECIST_single[[ref_name[input$show_details_single$index]]]
    # tmp$Plot = NA
    reactable( tmp,
               searchable = TRUE,
               paginationType = "jump",
               resizable = TRUE,
               showPageSizeOptions = TRUE,
               highlight = TRUE,
               bordered = TRUE,
               outlined = TRUE,
               striped = TRUE,
               defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                      sortNALast = TRUE,
                                      align = "center"
                                      # cell = function(value) format(value,digits = 2)
               ),
               # columns = list(
               #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
               #   `midp.exact`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `fisher.exact`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `chi.square`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `OR`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
               #     list(color = color, fontWeight = "bold")}
               #   )
               # ),
               # onClick = JS(
               #   "function(rowInfo, column){
               #                      if (column.id !== 'Plot') {return}
               #                      window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
               #                      if (window.Shiny) {
               #                                          Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
               #                                        }
               #           }"
               # ),
               defaultSortOrder = "asc",
               defaultSorted = list(chi.square = "asc", OR = "asc")

    )

  })
  output$RESPONSE_single = renderReactable({

    validate(need(!ref_name[input$show_details_single$index] %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                                             "dataset1","dataset13","dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21",
                                                             "dataset22","dataset23"),message = FALSE))
    tmp = ref_total_RESPONSE_single[[ref_name[input$show_details_single$index]]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `midp.exact`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   ),
                #   `fisher.exact`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   ),
                #   `chi.square`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   ),
                #   `OR`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(chi.square = "asc", OR = "asc")

    )
  })


  observeEvent(input$show_details_single$index,{

    ds = ref_name[input$show_details_single$index]

    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_explore_single_inf4")
      closeAlert(session,"warning_explore_single_sig4")
      createAlert(session, "warning_explore_single_inf", "warning_explore_single_inf4", title = "Warning",style = "danger",
                  content = paste("The dataset",ds,"didn't provide RNA-seq data",sep = " "), append = FALSE)
      createAlert(session, "warning_explore_single_sig", "warning_explore_single_sig4", title = "Warning",style = "danger",
                  content = paste("The dataset",ds,"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_explore_single_inf4")
      closeAlert(session,"warning_explore_single_sig4")
    }


  })
  output$immune_infiltrating_single = renderReactable({

    validate(need(ref_name[input$show_details_single$index] %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    tmp = ref_total_immune_infiltrating_single[[ref_name[input$show_details_single$index]]][[input$cell_type_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })
  output$immune_signature_single = renderReactable({

    validate(need(ref_name[input$show_details_single$index] %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    tmp = ref_total_immune_pathway_single[[ref_name[input$show_details_single$index]]][[input$signature_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )

  })
  ###################################Ref_datasets_Explore_pm######################################################################
  output$ref_datasets_detial_pm = renderReactable({
    datasets_overview$details = NA
    reactable(
      datasets_overview,
      rownames = F,
      searchable = TRUE,
      paginationType = "jump",
      resizable = TRUE,
      showPageSizeOptions = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      outlined = TRUE,
      striped = TRUE,
      defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                             sortNALast = TRUE,
                             align = "center"
                             # cell = function(value) format(value,digits = 2)
      ),
      columns = list(
        # Render a "show details" button in the last column of the table.
        # This button won't do anything by itself, but will trigger the custom
        # click action on the column.
        details = colDef(
          name = "",
          sortable = FALSE,
          cell = function() htmltools::tags$button("Show details")
        )
      ),
      onClick = JS("function(rowInfo, column) {
      if (column.id !== 'details') {
        return
      }
      if (window.Shiny) {
        Shiny.setInputValue('show_details_pm', { index: rowInfo.index + 1 }, { priority: 'event' })
        }
                     }"
      ),
      defaultSorted = list(Sequencing = "desc")
    )
  })
  observeEvent(input$show_details_pm$index,{

    ds = ref_name[input$show_details_pm$index]

    if(ds %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"]){
      closeAlert(session,"warning_explore_pm_OS2")
      closeAlert(session,"warning_explore_pm_PFS2")
      createAlert(session, "warning_explore_pm_OS", "warning_explore_pm_OS2", title = "Tips",
                  content = paste(ds,"consists of Panel sequencing data, which is not suitable for pathway mutation analysis."), append = FALSE)
      createAlert(session, "warning_explore_pm_PFS", "warning_explore_pm_PFS2", title = "Tips",
                  content = paste(ds,"consists of Panel sequencing data, which is not suitable for pathway mutation analysis."), append = FALSE)
    }else if(ds %in% c("dataset5","dataset7","dataset13","dataset16")){
      closeAlert(session,"warning_explore_pm_OS2")
      createAlert(session, "warning_explore_pm_OS", "warning_explore_pm_OS2", title = "Tips",
                  content = paste(ds,"didn't provide overall survival"), append = FALSE)
      closeAlert(session,"warning_explore_pm_PFS2")
    }else if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                       "dataset1","dataset8","dataset10","dataset11")){
      closeAlert(session,"warning_explore_pm_PFS2")
      createAlert(session, "warning_explore_pm_PFS", "warning_explore_pm_PFS2", title = "Tips",
                  content = paste(ds,"didn't provide progression-free survival"), append = FALSE)
      closeAlert(session,"warning_explore_pm_OS2")
    }else{
      closeAlert(session,"warning_explore_pm_OS2")
      closeAlert(session,"warning_explore_pm_PFS2")
    }

  })
  output$OS_total_pm = renderReactable({

    validate(need(!ref_name[input$show_details_pm$index] %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"],message = FALSE))
    validate(need(!ref_name[input$show_details_pm$index] %in% c("dataset5","dataset7","dataset13","dataset16"),message = FALSE))
    tmp = ref_total_OS_pm[[ref_name[input$show_details_pm$index]]]
    # tmp$Plot = NA
    reactable( tmp,
               searchable = TRUE,
               paginationType = "jump",
               resizable = TRUE,
               showPageSizeOptions = TRUE,
               #             onClick = "expand",
               highlight = TRUE,
               bordered = TRUE,
               outlined = TRUE,
               striped = TRUE,
               defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                      sortNALast = TRUE,
                                      align = "center"
                                      # cell = function(value) format(value,digits = 2)
                                      ),
               # columns = list(
               #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
               #   `Log_rank_test(OS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `Wald_test(OS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `HR(OS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
               #     list(color = color, fontWeight = "bold")}
               #   )
               # ),
               # onClick = JS(
               #   "function(rowInfo, column){
               #                      if (column.id !== 'Plot') {return}
               #                      window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
               #                      if (window.Shiny) {
               #                                          Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
               #                                        }
               #           }"
               # ),
               defaultSortOrder = "asc",
               defaultSorted = list(`Log_rank_test(OS)` = "asc", `HR(OS)` = "asc")

    )
  })
  output$PFS_total_pm = renderReactable({

    validate(need(!ref_name[input$show_details_pm$index] %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"],message = FALSE))
    validate(need(!ref_name[input$show_details_pm$index] %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                                                               "dataset1","dataset8","dataset10","dataset11"),message = FALSE))

    tmp = ref_total_PFS_pm[[ref_name[input$show_details_pm$index]]]
    # tmp$Plot = NA
    reactable( tmp,
               searchable = TRUE,
               paginationType = "jump",
               resizable = TRUE,
               showPageSizeOptions = TRUE,
               #             onClick = "expand",
               highlight = TRUE,
               bordered = TRUE,
               outlined = TRUE,
               striped = TRUE,
               defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                      sortNALast = TRUE,
                                      align = "center"
                                      # cell = function(value) format(value,digits = 2)
                                      ),
               # columns = list(
               #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
               #   `Log_rank_test(PFS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `Wald_test(PFS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `HR(PFS)`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
               #     list(color = color, fontWeight = "bold")}
               #   )
               # ),
               # onClick = JS(
               #   "function(rowInfo, column){
               #                      if (column.id !== 'Plot') {return}
               #                      window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
               #                      if (window.Shiny) {
               #                                          Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
               #                                        }
               #           }"
               # ),
               defaultSortOrder = "asc",
               defaultSorted = list(`Log_rank_test(PFS)` = "asc", `HR(PFS)` = "asc")

    )
  })

  observeEvent(input$show_details_pm$index,{

    ds = ref_name[input$show_details_pm$index]

    if(ds %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"]){
      closeAlert(session,"warning_explore_pm_recist2")
      closeAlert(session,"warning_explore_pm_res2")
      createAlert(session, "warning_explore_pm_recist", "warning_explore_pm_recist2", title = "Tips",
                  content = paste(ds,"consists of Panel sequencing data, which is not suitable for pathway mutation analysis."), append = FALSE)
      createAlert(session, "warning_explore_pm_res", "warning_explore_pm_res2", title = "Tips",
                  content = paste(ds,"consists of Panel sequencing data, which is not suitable for pathway mutation analysis."), append = FALSE)
    }else if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                 "dataset1","dataset13")){
      closeAlert(session,"warning_explore_pm_recist2")
      closeAlert(session,"warning_explore_pm_res2")
      createAlert(session, "warning_explore_pm_recist", "warning_explore_pm_recist2", title = "Warning",style = "danger",
                  content = paste(ds,"didn't provide data about RECIST"), append = FALSE)
      createAlert(session, "warning_explore_pm_res", "warning_explore_pm_res2", title = "Warning",style = "danger",
                  content = paste(ds,"didn't provide data about DCB"), append = FALSE)

    }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
      closeAlert(session,"warning_explore_pm_recist2")
      createAlert(session, "warning_explore_pm_recist", "warning_explore_pm_recist2", title = "Tips",
                  content = paste(ds,"didn't provide data about RECIST"), append = FALSE)
      closeAlert(session,"warning_explore_pm_res2")
    }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
      closeAlert(session,"warning_explore_pm_res2")
      createAlert(session, "warning_explore_pm_res", "warning_explore_pm_res2", title = "Tips",
                  content = paste(ds,"didn't provide data about DCB"), append = FALSE)
      closeAlert(session,"warning_explore_pm_recist2")
    }else{
      closeAlert(session,"warning_explore_pm_recist2")
      closeAlert(session,"warning_explore_pm_res2")
    }



  })
  output$RECIST_pm = renderReactable({

    validate(need(!ref_name[input$show_details_pm$index] %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"],message = FALSE))
    validate(need(!ref_name[input$show_details_pm$index] %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                                                               "dataset1","dataset13","dataset5","dataset6","dataset9","dataset10","dataset15","dataset19"),message = FALSE))

    tmp = ref_total_RECIST_pm[[ref_name[input$show_details_pm$index]]]
    # tmp$Plot = NA
    reactable( tmp,
               searchable = TRUE,
               paginationType = "jump",
               resizable = TRUE,
               showPageSizeOptions = TRUE,
               highlight = TRUE,
               bordered = TRUE,
               outlined = TRUE,
               striped = TRUE,
               defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                      sortNALast = TRUE,
                                      align = "center"
                                      # cell = function(value) format(value,digits = 2)
                                      ),
               # columns = list(
               #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
               #   `midp.exact`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `fisher.exact`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `chi.square`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
               #     list(color = color, fontWeight = "bold")}
               #   ),
               #   `OR`=colDef(style = function(value){
               #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
               #     list(color = color, fontWeight = "bold")}
               #   )
               # ),
               # onClick = JS(
               #   "function(rowInfo, column){
               #                      if (column.id !== 'Plot') {return}
               #                      window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
               #                      if (window.Shiny) {
               #                                          Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
               #                                        }
               #           }"
               # ),
               defaultSortOrder = "asc",
               defaultSorted = list(chi.square = "asc", OR = "asc")

    )
  })
  output$RESPONSE_pm = renderReactable({

    validate(need(!ref_name[input$show_details_pm$index] %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"],message = FALSE))
    validate(need(!ref_name[input$show_details_pm$index] %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                                                               "dataset1","dataset13","dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21",
                                                                               "dataset22","dataset23"),message = FALSE))

    tmp = ref_total_RESPONSE_pm[[ref_name[input$show_details_pm$index]]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `midp.exact`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   ),
                #   `fisher.exact`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   ),
                #   `chi.square`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   ),
                #   `OR`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<1){color <- "blue"}else{color <- "red"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(chi.square = "asc", OR = "asc")

    )
  })


  observeEvent(input$show_details_pm$index,{

    ds = ref_name[input$show_details_pm$index]

    if(ds %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"]){

      closeAlert(session,"warning_explore_pm_inf4")
      closeAlert(session,"warning_explore_pm_sig4")
      createAlert(session, "warning_explore_pm_inf", "warning_explore_pm_inf4", title = "Tips",
                  content = paste(ds,"consists of Panel sequencing data, which is not suitable for pathway mutation analysis."), append = FALSE)
      createAlert(session, "warning_explore_pm_sig", "warning_explore_pm_sig4", title = "Tips",
                  content = paste(ds,"consists of Panel sequencing data, which is not suitable for pathway mutation analysis."), append = FALSE)

    }else if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_explore_pm_inf4")
      closeAlert(session,"warning_explore_pm_sig4")
      createAlert(session, "warning_explore_pm_inf", "warning_explore_pm_inf4", title = "Warning",style = "danger",
                  content = paste("The dataset",ds,"didn't provide RNA-seq data",sep = " "), append = FALSE)
      createAlert(session, "warning_explore_pm_sig", "warning_explore_pm_sig4", title = "Warning",style = "danger",
                  content = paste("The dataset",ds,"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_explore_pm_inf4")
      closeAlert(session,"warning_explore_pm_sig4")
    }


  })
  output$immune_infiltrating_pm = renderReactable({

    validate(need(!ref_name[input$show_details_pm$index] %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"],message = FALSE))
    validate(need(ref_name[input$show_details_pm$index] %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    tmp = ref_total_immune_infiltrating_pm[[ref_name[input$show_details_pm$index]]][[input$cell_type_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })
  output$immune_signature_pm = renderReactable({

    validate(need(!ref_name[input$show_details_pm$index] %in% rownames(datasets_overview)[datasets_overview$Sequencing == "Panel"],message = FALSE))
    validate(need(ref_name[input$show_details_pm$index] %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))


    tmp = ref_total_immune_pathway_pm[[ref_name[input$ref_datasets_detial_pm_rows_selected]]][[input$signature_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  # observeEvent(input$ref_datasets_detial_rows_selected,{message(datasets_overview[input$ref_datasets_detial_rows_selected,])})
  ###################################TCGA_Explore_single######################################################################
  output$TCGA_detial_single = renderReactable({
    TCGA_overview$details = NA
    reactable(
      TCGA_overview,
      rownames = F,
      searchable = TRUE,
      paginationType = "jump",
      resizable = TRUE,
      showPageSizeOptions = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      outlined = TRUE,
      striped = TRUE,
      defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                             sortNALast = TRUE,
                             align = "center"
                             # cell = function(value) format(value,digits = 2)
      ),
      columns = list(
        # Render a "show details" button in the last column of the table.
        # This button won't do anything by itself, but will trigger the custom
        # click action on the column.
        details = colDef(
          name = "",
          sortable = FALSE,
          cell = function() htmltools::tags$button("Show details")
        )
      ),
      onClick = JS("function(rowInfo, column) {
      if (column.id !== 'details') {
        return
      }
      if (window.Shiny) {
        Shiny.setInputValue('show_details_tcga', { index: rowInfo.index + 1 }, { priority: 'event' })
        }
                     }"
      )
    )
  })

  output$TCGA_immune_infiltrating_single = renderReactable({
    validate(need(TCGA_name[input$show_details_tcga$index] %in% names(TCGA),message = FALSE))
    tmp = TCGA_total_immune_infiltration_single[[TCGA_name[input$show_details_tcga$index]]][[input$TCGA_cell_type_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$TCGA_immune_signature_single = renderReactable({
    validate(need(TCGA_name[input$show_details_tcga$index] %in% names(TCGA),message = FALSE))
    tmp = TCGA_total_immune_pathway_single[[TCGA_name[input$show_details_tcga$index]]][[input$TCGA_signature_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })


  ###################################TCGA_Explore_pm######################################################################
  output$TCGA_detial_pm = renderReactable({
    TCGA_overview$details = NA
    reactable(
      TCGA_overview,
      rownames = F,
      searchable = TRUE,
      paginationType = "jump",
      resizable = TRUE,
      showPageSizeOptions = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      outlined = TRUE,
      striped = TRUE,
      defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                             sortNALast = TRUE,
                             align = "center"
                             # cell = function(value) format(value,digits = 2)
      ),
      columns = list(
        # Render a "show details" button in the last column of the table.
        # This button won't do anything by itself, but will trigger the custom
        # click action on the column.
        details = colDef(
          name = "",
          sortable = FALSE,
          cell = function() htmltools::tags$button("Show details")
        )
      ),
      onClick = JS("function(rowInfo, column) {
      if (column.id !== 'details') {
        return
      }
      if (window.Shiny) {
        Shiny.setInputValue('show_details_tcga_pm', { index: rowInfo.index + 1 }, { priority: 'event' })
        }
                     }"
      )
    )
  })
  output$TCGA_immune_infiltrating_pm = renderReactable({
    validate(need(TCGA_name[input$show_details_tcga_pm$index] %in% names(TCGA),message = FALSE))
    tmp = TCGA_total_immune_infiltration_pm[[TCGA_name[input$show_details_tcga_pm$index]]][[input$TCGA_cell_type_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$TCGA_immune_signature_pm = renderReactable({
    validate(need(TCGA_name[input$show_details_tcga_pm$index] %in% names(TCGA),message = FALSE))
    tmp = TCGA_total_immune_pathway_pm[[TCGA_name[input$show_details_tcga_pm$index]]][[input$TCGA_signature_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  ###################################CPTAC_Explore_single######################################################################
  output$CPTAC_detial_single = renderReactable({
    CPTAC_overview$details = NA
    reactable(
      CPTAC_overview,
      rownames = F,
      searchable = TRUE,
      paginationType = "jump",
      resizable = TRUE,
      showPageSizeOptions = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      outlined = TRUE,
      striped = TRUE,
      defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                             sortNALast = TRUE,
                             align = "center"
                             # cell = function(value) format(value,digits = 2)
      ),
      columns = list(
        # Render a "show details" button in the last column of the table.
        # This button won't do anything by itself, but will trigger the custom
        # click action on the column.
        details = colDef(
          name = "",
          sortable = FALSE,
          cell = function() htmltools::tags$button("Show details")
        )
      ),
      onClick = JS("function(rowInfo, column) {
      if (column.id !== 'details') {
        return
      }
      if (window.Shiny) {
        Shiny.setInputValue('show_details_cptac', { index: rowInfo.index + 1 }, { priority: 'event' })
        }
                     }"
      )
    )
  })

  output$CPTAC_immune_infiltrating_rna_single = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_infiltration_rna_single[[CPTAC_name[input$show_details_cptac$index]]][[input$CPTAC_cell_type_rna_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$CPTAC_immune_signature_rna_single = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_pathway_rna_single[[CPTAC_name[input$show_details_cptac$index]]][[input$CPTAC_signature_rna_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$CPTAC_immune_infiltrating_protein_single = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_infiltration_protein_single[[CPTAC_name[input$show_details_cptac$index]]][[input$CPTAC_cell_type_protein_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$CPTAC_immune_signature_protein_single = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_pathway_protein_single[[CPTAC_name[input$show_details_cptac$index]]][[input$CPTAC_signature_protein_single]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  ###################################CPTAC_Explore_pm######################################################################
  output$CPTAC_detial_pm = renderReactable({
    CPTAC_overview$details = NA
    reactable(
      CPTAC_overview,
      rownames = F,
      searchable = TRUE,
      paginationType = "jump",
      resizable = TRUE,
      showPageSizeOptions = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      outlined = TRUE,
      striped = TRUE,
      defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                             sortNALast = TRUE,
                             align = "center"
                             # cell = function(value) format(value,digits = 2)
      ),
      columns = list(
        # Render a "show details" button in the last column of the table.
        # This button won't do anything by itself, but will trigger the custom
        # click action on the column.
        details = colDef(
          name = "",
          sortable = FALSE,
          cell = function() htmltools::tags$button("Show details")
        )
      ),
      onClick = JS("function(rowInfo, column) {
      if (column.id !== 'details') {
        return
      }
      if (window.Shiny) {
        Shiny.setInputValue('show_details_cptac_pm', { index: rowInfo.index + 1 }, { priority: 'event' })
        }
                     }"
      )
    )
  })

  output$CPTAC_immune_infiltrating_rna_pm = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac_pm$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_infiltration_rna_pm[[CPTAC_name[input$show_details_cptac_pm$index]]][[input$CPTAC_cell_type_rna_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$CPTAC_immune_signature_rna_pm = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac_pm$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_pathway_rna_pm[[CPTAC_name[input$show_details_cptac_pm$index]]][[input$CPTAC_signature_rna_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$CPTAC_immune_infiltrating_protein_pm = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac_pm$index] %in% names(CPTAC),message = FALSE))
    message(input$CPTAC_cell_type_protein_pm)
    message(CPTAC_name[input$show_details_cptac_pm$index])
    tmp = CPTAC_total_immune_infiltration_protein_pm[[CPTAC_name[input$show_details_cptac_pm$index]]][[input$CPTAC_cell_type_protein_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  output$CPTAC_immune_signature_protein_pm = renderReactable({
    validate(need(CPTAC_name[input$show_details_cptac_pm$index] %in% names(CPTAC),message = FALSE))
    tmp = CPTAC_total_immune_pathway_protein_pm[[CPTAC_name[input$show_details_cptac_pm$index]]][[input$CPTAC_signature_protein_pm]]
    # tmp$Plot = NA
    reactable(  tmp,
                searchable = TRUE,
                paginationType = "jump",
                resizable = TRUE,
                showPageSizeOptions = TRUE,
                highlight = TRUE,
                bordered = TRUE,
                outlined = TRUE,
                striped = TRUE,
                defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                       sortNALast = TRUE,
                                       align = "center"
                                       # cell = function(value) format(value,digits = 2)
                                       ),
                # columns = list(
                #   Plot=colDef(cell = function() htmltools::tags$button("Show details")),
                #   `P_value`=colDef(style = function(value){
                #     if(is.na(value)){color <- "black"}else if(value<0.05){color <- "#e00000"}else{color <- "black"}
                #     list(color = color, fontWeight = "bold")}
                #   )
                #
                # ),
                # onClick = JS(
                #   "function(rowInfo, column){
                #                         if (column.id !== 'Plot') {return}
                #                         window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                #                         if (window.Shiny) {
                #                                             Shiny.setInputValue('show_details', { index: {index:rowInfo.index + 1,values:rowInfo.values} }, { priority: 'event' })
                #                                           }
                #              }"
                # ),
                defaultSortOrder = "asc",
                defaultSorted = list(P_value = "asc")

    )
  })

  ####################################Search single gene###############################################
  search_gene = eventReactive(input$search_single,{input$search_gene})

  output$gene_detial = renderUI({HTML(Gene_detial_information(search_gene()))})
  output$OS_single_forest = renderPlot({OS_single_forest(search_gene())})
  output$PFS_single_forest = renderPlot({PFS_single_forest(search_gene())})
  output$RECIST_single_forest = renderPlot({RECIST_single_forest(search_gene())})
  output$RESPONSE_single_forest = renderPlot({RESPONSE_single_forest(search_gene())})
  output$ref_single_Immune_Infiltration_heatmap = renderPlot({ref_single_Immune_Infiltration_heatmap(search_gene())})
  output$ref_single_Immune_pathway_heatmap = renderPlot({ref_single_Immune_pathway_heatmap(search_gene())})
  output$TCGA_single_Immune_Infiltration_heatmap = renderPlot({TCGA_single_Immune_Infiltration_heatmap(search_gene())})
  output$TCGA_single_Immune_pathway_heatmap = renderPlot({TCGA_single_Immune_pathway_heatmap(search_gene())})

  output$CPTAC_rna_single_Immune_Infiltration_heatmap = renderPlot({CPTAC_rna_single_Immune_Infiltration_heatmap(search_gene())})
  output$CPTAC_protein_single_Immune_Infiltration_heatmap = renderPlot({CPTAC_protein_single_Immune_Infiltration_heatmap(search_gene())})
  output$CPTAC_rna_single_Immune_pathway_heatmap = renderPlot({CPTAC_rna_single_Immune_pathway_heatmap(search_gene())})
  output$CPTAC_protein_single_Immune_pathway_heatmap = renderPlot({CPTAC_protein_single_Immune_pathway_heatmap(search_gene())})
  ####################################Search pathway###############################################
  observeEvent(input$search_pm_kw,{
    target_pathway = names(pathway_list)[grepl(input$search_pathway_kw,names(pathway_list))]
    updatePickerInput(session,inputId = "search_pathway",choices = target_pathway,selected = target_pathway[1])
  })

  search_pathway = eventReactive(input$search_pm,{input$search_pathway})

  output$OS_pm_forest = renderPlot({OS_pm_forest(search_pathway())})
  output$PFS_pm_forest = renderPlot({PFS_pm_forest(search_pathway())})
  output$RECIST_pm_forest = renderPlot({RECIST_pm_forest(search_pathway())})
  output$RESPONSE_pm_forest = renderPlot({RESPONSE_pm_forest(search_pathway())})
  output$ref_pm_Immune_Infiltration_heatmap = renderPlot({ref_pm_Immune_Infiltration_heatmap(search_pathway())})
  output$ref_pm_Immune_pathway_heatmap = renderPlot({ref_pm_Immune_pathway_heatmap(search_pathway())})
  output$TCGA_pm_Immune_Infiltration_heatmap = renderPlot({TCGA_pm_Immune_Infiltration_heatmap(search_pathway())})
  output$TCGA_pm_Immune_pathway_heatmap = renderPlot({TCGA_pm_Immune_pathway_heatmap(search_pathway())})

  output$CPTAC_rna_pm_Immune_Infiltration_heatmap = renderPlot({CPTAC_rna_pm_Immune_Infiltration_heatmap(search_pathway())})
  output$CPTAC_protein_pm_Immune_Infiltration_heatmap = renderPlot({CPTAC_protein_pm_Immune_Infiltration_heatmap(search_pathway())})
  output$CPTAC_rna_pm_Immune_pathway_heatmap = renderPlot({CPTAC_rna_pm_Immune_pathway_heatmap(search_pathway())})
  output$CPTAC_protein_pm_Immune_pathway_heatmap = renderPlot({CPTAC_protein_pm_Immune_pathway_heatmap(search_pathway())})
}

shinyApp(ui = ui, server = server)


