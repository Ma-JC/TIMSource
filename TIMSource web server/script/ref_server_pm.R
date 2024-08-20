Ref_datasets_server_pm <- function(input,output,session,ref_pm_para,datasets,datasets_mu,min.pct_ref_pm,FC_ref_pm,pvalue_ref_pm,useid){
  
  useid <- paste(useid,"Refpm",sep = "_")
  
  width = reactive({if(ref_pm_para()$datasource_pm %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset4","dataset6","dataset9","dataset12","dataset14","dataset15","dataset17","dataset18","dataset19","dataset20","dataset21","dataset22","dataset23","dataset24")){900}else{500}})
  width2 = reactive({if(ref_pm_para()$datasource_pm %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset7","dataset12","dataset14","dataset16","dataset24")){1200}else{600}})
  width3 = reactive({if(ref_pm_para()$datasource_pm %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9","dataset3.1","dataset3.2","dataset3.3","dataset3",
                                               "dataset1","dataset2","dataset4","dataset5","dataset6","dataset7","dataset8","dataset9","dataset10","dataset12","dataset14","dataset15","dataset16","dataset17","dataset20","dataset21")){1600}else{1000}})
  
  
  ref_cohort_pm = reactive({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    Mut_type_ref_pm = ref_pm_para()$Mut_type_ref_pm
    Wild_type_ref_pm = ref_pm_para()$Wild_type_ref_pm
    Therapy_type_ref_pm = ref_pm_para()$Therapy_type_ref_pm
    
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    ref_cohort_cal_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],Mut_type = Mut_type_ref_pm,Wild_type = Wild_type_ref_pm,therapy_type=Therapy_type_ref_pm)
  }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_pm_para()$Mut_type_ref_pm,ref_pm_para()$Wild_type_ref_pm,ref_pm_para()$Therapy_type_ref_pm)
  
  #################################################### survival #################################
  output$sur_pm_uidown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(!ds %in% c(),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_sur_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_sur_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sur_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sur_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_sur_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$ref_sur_pm_file,{
    if(input$ref_sur_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_sur_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_sur_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_sur_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_sur_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_sur_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_sur_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    if(ds %in% c("dataset5","dataset7","dataset13","dataset16")){
      closeAlert(session,"warning_id_pm")
      createAlert(session, "warning_pm", "warning_id_pm", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide overall survival"), append = FALSE)
    }else if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                       "dataset1","dataset8","dataset10","dataset11")){
      closeAlert(session,"warning_id_pm")
      createAlert(session, "warning_pm", "warning_id_pm", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide progression-free survival"), append = FALSE)
    }else{

      closeAlert(session,"warning_id_pm")
    }
    
    if(!any(pathway_list[[gene]] %in% colnames(datasets[[ds]]))){
      closeAlert(session,"warning_id_pm")
      createAlert(session, "warning_pm", "warning_id_pm", title = "Warning",style = "danger",
                  content = paste("The genes of",gene,"does not exist in",dataset_name2[[ds]],sep = " "), append = FALSE)
    }
    
  })
  
  observe({
    gene = ref_pm_para()$gene_pm
    ref_cohort_pm = ref_cohort_pm()
    
    if(length(ref_cohort_pm$mut) <3 | length(ref_cohort_pm$wt) < 3){
      closeAlert(session,"warning_id_pm")
      createAlert(session, "warning_pm", "warning_id_pm", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }
    
  })
  
  output$Survival_pm = renderPlot({
    
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    width = width()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    #req(iv$is_valid())
    if(ds %in%  c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset4","dataset6","dataset9","dataset12","dataset14",
                  "dataset15","dataset17","dataset18","dataset19","dataset20","dataset21","dataset22","dataset23","dataset24")){
      p1 = os_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
      p2 = pfs_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
      layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
      tmp = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)
    }else if(ds %in%  c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                        "dataset1","dataset8","dataset10","dataset11")){
      tmp = os_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    }else if(ds %in%  c("dataset5","dataset7","dataset13","dataset16")){
      tmp = pfs_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    }
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)
    
  },width = width,height = 600)
  
  output$ref_sur_pm_down = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      
      if(input$ref_sur_pm_file == 'pdf'){
        paste0(gene,"_",input$dataset_pm,".","sur", '.',"pdf")
      }else if(input$ref_sur_pm_file == 'png'){
        paste0(gene,"_",input$dataset_pm,".","sur", '.',"png")
      }else if(input$ref_sur_pm_file == 'jpeg'){
        paste0(gene,"_",input$dataset_pm,".","sur", '.',"jpeg")
      }else{
        paste0(gene,"_",input$dataset_pm,".","sur", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("ref_sur_pm_down")
      gene = ref_pm_para()$gene_pm
      ds = ref_pm_para()$datasource_pm
      ref_cohort_pm = ref_cohort_pm()
      
      if(input$ref_sur_pm_res <= 300){
        r = input$ref_sur_pm_res
      }else{
        r = 300
      }
      if(input$ref_sur_pm_file == "pdf" & input$ref_sur_pm_height <= 50){
        h = input$ref_sur_pm_height
      }else if(input$ref_sur_pm_file == "pdf" & input$ref_sur_pm_height > 50){
        h = 50
      }else if(input$ref_sur_pm_file != "pdf" & input$ref_sur_pm_height <= 3000){
        h = input$ref_sur_pm_height
      }else{
        h = 3000
      }
      if(input$ref_sur_pm_file == "pdf" & input$ref_sur_pm_width <= 50){
        w = input$ref_sur_pm_width
      }else if(input$ref_sur_pm_file == "pdf" & input$ref_sur_pm_width > 50){
        w = 50
      }else if(input$ref_sur_pm_file != "pdf" & input$ref_sur_pm_width <= 3000){
        w = input$ref_sur_pm_width
      }else{
        w = 3000
      }
      
      if(input$ref_sur_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_sur_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_sur_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      if(ds %in%  c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset4","dataset6","dataset9","dataset12","dataset14",
                    "dataset15","dataset17","dataset18","dataset19","dataset20","dataset21","dataset22","dataset23","dataset24")){
        p1 = os_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
        p2 = pfs_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
        layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
        print((p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout))
      }else if(ds %in%  c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                          "dataset1","dataset8","dataset10","dataset11")){
        print(os_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      }else if(ds %in%  c("dataset5","dataset7","dataset13","dataset16")){
        print(pfs_survival_pm(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      }
      dev.off()
      shinyjs::enable("ref_sur_pm_down")
    }
  )
  
  # output$warning_pm = renderText({
  #   gene = ref_pm_para()$gene_pm
  #   ds = ref_pm_para()$datasource_pm
  #   if(ds %in% c("dataset5","dataset7","dataset13","dataset16")){
  #     paste(ds,"didn't provide overall survival")
  #   }else if(ds %in% c("dataset1","dataset8","dataset10","dataset11")){
  #     paste(ds,"didn't provide progression-free survival")
  #   }
  #   
  # })
  
  ######################################################  Response  ################################################
  
  output$res_pm_uidown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(!ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                             "dataset1","dataset13"),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_res_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_res_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_res_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_res_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_res_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$ref_res_pm_file,{
    if(input$ref_res_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_res_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_res_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_res_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_res_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_res_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_res_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    
    if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                 "dataset1","dataset13")){
      closeAlert(session,"warning_id2_pm")
      createAlert(session, "warning2_pm", "warning_id2_pm", title = "Warning",style = "danger",
                  content = paste(dataset_name2[[ds]],"did not provide any data on the efficacy of immunotherapy"), append = FALSE)
      
    }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
      closeAlert(session,"warning_id2_pm")
      createAlert(session, "warning2_pm", "warning_id2_pm", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide data about RECIST"), append = FALSE)
      
    }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
      closeAlert(session,"warning_id2_pm")
      createAlert(session, "warning2_pm", "warning_id2_pm", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide data about DCB"), append = FALSE)
    }else{
      closeAlert(session,"warning_id2_pm")
    }
    
    if(!any(pathway_list[[gene]] %in% colnames(datasets[[ds]]))){
      closeAlert(session,"warning_id2_pm")
      createAlert(session, "warning2_pm", "warning_id2_pm", title = "Warning",style = "danger",
                  content = paste("The genes of",gene,"does not exist in",dataset_name2[[ds]],sep = " "), append = FALSE)
    }
    
  })
  
  observe({
    gene = ref_pm_para()$gene_pm
    ref_cohort_pm = ref_cohort_pm()
    
    if(length(ref_cohort_pm$mut) <3 | length(ref_cohort_pm$wt) < 3){
      closeAlert(session,"warning_id2_pm")
      createAlert(session, "warning2_pm", "warning_id2_pm", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }
    
  })
  
  output$response_pm = renderPlot({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    width2 = width2()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    validate(need(!ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                             "dataset1","dataset13"),message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    if(ds %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset7","dataset12","dataset14","dataset16","dataset24")){
      p1 = CRPR_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
      p2 = DCB_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
      tmp = p1 + p2
    }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
      tmp = DCB_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
      tmp = CRPR_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    }
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)

  },width = width2,height = 600)
  
  output$ref_res_pm_down = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      
      if(input$ref_res_pm_file == 'pdf'){
        paste0(gene,"_",input$dataset_pm,".","res", '.',"pdf")
      }else if(input$ref_res_pm_file == 'png'){
        paste0(gene,"_",input$dataset_pm,".","res", '.',"png")
      }else if(input$ref_res_pm_file == 'jpeg'){
        paste0(gene,"_",input$dataset_pm,".","res", '.',"jpeg")
      }else{
        paste0(gene,"_",input$dataset_pm,".","res", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("ref_res_pm_down")
      gene = ref_pm_para()$gene_pm
      ds = ref_pm_para()$datasource_pm
      ref_cohort_pm = ref_cohort_pm()
      
      if(input$ref_res_pm_res <= 300){
        r = input$ref_res_pm_res
      }else{
        r = 300
      }
      if(input$ref_res_pm_file == "pdf" & input$ref_res_pm_height <= 50){
        h = input$ref_res_pm_height
      }else if(input$ref_res_pm_file == "pdf" & input$ref_res_pm_height > 50){
        h = 50
      }else if(input$ref_res_pm_file != "pdf" & input$ref_res_pm_height <= 3000){
        h = input$ref_res_pm_height
      }else{
        h = 3000
      }
      if(input$ref_res_pm_file == "pdf" & input$ref_res_pm_width <= 50){
        w = input$ref_res_pm_width
      }else if(input$ref_res_pm_file == "pdf" & input$ref_res_pm_width > 50){
        w = 50
      }else if(input$ref_res_pm_file != "pdf" & input$ref_res_pm_width <= 3000){
        w = input$ref_res_pm_width
      }else{
        w = 3000
      }
      
      if(input$ref_res_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_res_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_res_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      if(ds %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset7","dataset12","dataset14","dataset16","dataset24")){
        p1 = CRPR_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
        p2 = DCB_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
        print(p1 + p2)
      }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
        print(DCB_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
        print(CRPR_pm(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      }
      dev.off()
      shinyjs::enable("ref_res_pm_down")
    }
  )
  
  # output$warning2_pm = renderText({
  #   gene = ref_pm_para()$gene_pm
  #   ds = ref_pm_para()$datasource_pm
  #   validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),paste("The genes of",gene,"does not exist in",ds,sep = " ")))
  #   if(ds %in% c("dataset1","dataset13")){
  #     "Data set 1 did not provide any data on the efficacy of immunotherapy"
  #   }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
  #     paste(ds,"didn't provide data about RECIST")
  #   }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20")){
  #     paste(ds,"didn't provide data about DCB")
  #   }
  #   
  # })
  
  ######################################################### Mutation #######################################
  
  output$mut_pm_uidown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(!ds %in% c(),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_mut_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_mut_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_mut_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_mut_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_mut_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$ref_mut_pm_file,{
    if(input$ref_mut_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_mut_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_mut_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_mut_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_mut_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_mut_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_mut_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_pm_uitabledown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(!ds %in% c(),message = FALSE))
    downloadButton('ref_mut_pm_tabdown',label = 'Download Table')
  })
  
  
  observe({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm

    if(ds %in% c("dataset18","dataset19","dataset22","dataset23")){
      closeAlert(session,"warning_id3_pm")
      createAlert(session, "warning3_pm", "warning_id3_pm", title = "Tips",
                  content = paste(dataset_name2[[ds]],"did not provide TMB data"), append = FALSE)
      
    }else if(ds %in% c("dataset11","dataset13")){
      closeAlert(session,"warning_id3_pm")
      createAlert(session, "warning3_pm", "warning_id3_pm", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide detial information about gene mutation"), append = FALSE)
      
    }else{
      closeAlert(session,"warning_id3_pm")
    }
    
    if(!any(pathway_list[[gene]] %in% colnames(datasets[[ds]]))){
      closeAlert(session,"warning_id3_pm")
      createAlert(session, "warning3_pm", "warning_id3_pm", title = "Warning",style = "danger",
                  content = paste("The genes of",gene,"does not exist in",dataset_name2[[ds]],sep = " "), append = FALSE)
    }
    
  })
  
  
  observe({
    gene = ref_pm_para()$gene_pm
    ref_cohort_pm = ref_cohort_pm()
    
    if(length(ref_cohort_pm$mut) <3 | length(ref_cohort_pm$wt) < 3){
      closeAlert(session,"warning_id3_pm")
      createAlert(session, "warning3_pm", "warning_id3_pm", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }

  })
  
  
  output$TMB_pm = renderPlot({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    width3 = width3()
    outline = input$outline_pm
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    #validate(need(! ds %in% c("dataset10"),paste("The dataset",ds,"didn't provide TMB or detial mutation data ",sep = " ")))
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9","dataset3.1","dataset3.2","dataset3.3","dataset3",
                 "dataset1","dataset2","dataset4","dataset5","dataset6","dataset7","dataset8","dataset9","dataset10","dataset12","dataset14","dataset15","dataset16","dataset17","dataset20","dataset21")){
      
      if(outline){
        p1 = TMB_pm(datasets[[ds]],gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
      }else{p1 = TMB_pm(datasets[[ds]],gene,remove_outlier = FALSE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)}
      
      p2 = MutationType_pm(dataset_mu = datasets_mu[[ds]],gene = gene)
      tmp = p1 + p2
    }else if(ds %in% c("dataset18","dataset19","dataset22","dataset23")){
      tmp = MutationType_pm(dataset_mu = datasets_mu[[ds]],gene = gene)
    }else if(ds %in% c("dataset11","dataset13","dataset24")){
      tmp = TMB_pm(datasets[[ds]],gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    }
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)
  },width = width3,height = 600)
  
  output$ref_mut_pm_down = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      
      if(input$ref_mut_pm_file == 'pdf'){
        paste0(gene,"_",input$dataset_pm,".","mut", '.',"pdf")
      }else if(input$ref_mut_pm_file == 'png'){
        paste0(gene,"_",input$dataset_pm,".","mut", '.',"png")
      }else if(input$ref_mut_pm_file == 'jpeg'){
        paste0(gene,"_",input$dataset_pm,".","mut", '.',"jpeg")
      }else{
        paste0(gene,"_",input$dataset_pm,".","mut", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("ref_mut_pm_down")
      gene = ref_pm_para()$gene_pm
      ds = ref_pm_para()$datasource_pm
      ref_cohort_pm = ref_cohort_pm()
      
      if(input$ref_mut_pm_res <= 300){
        r = input$ref_mut_pm_res
      }else{
        r = 300
      }
      if(input$ref_mut_pm_file == "pdf" & input$ref_mut_pm_height <= 50){
        h = input$ref_mut_pm_height
      }else if(input$ref_mut_pm_file == "pdf" & input$ref_mut_pm_height > 50){
        h = 50
      }else if(input$ref_mut_pm_file != "pdf" & input$ref_mut_pm_height <= 3000){
        h = input$ref_mut_pm_height
      }else{
        h = 3000
      }
      if(input$ref_mut_pm_file == "pdf" & input$ref_mut_pm_width <= 50){
        w = input$ref_mut_pm_width
      }else if(input$ref_mut_pm_file == "pdf" & input$ref_mut_pm_width > 50){
        w = 50
      }else if(input$ref_mut_pm_file != "pdf" & input$ref_mut_pm_width <= 3000){
        w = input$ref_mut_pm_width
      }else{
        w = 3000
      }
      
      if(input$ref_mut_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_mut_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_mut_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9","dataset3.1","dataset3.2","dataset3.3","dataset3",
                   "dataset1","dataset2","dataset4","dataset5","dataset6","dataset7","dataset8","dataset9","dataset10","dataset12","dataset14","dataset15","dataset16","dataset17","dataset20","dataset21")){
        
        if(input$outline_pm){
          p1 = TMB_pm(datasets[[ds]],gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
        }else{p1 = TMB_pm(datasets[[ds]],gene,remove_outlier = FALSE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)}
        
        p2 = MutationType_pm(dataset_mu = datasets_mu[[ds]],gene = gene)
        print(p1 + p2)
      }else if(ds %in% c("dataset18","dataset19","dataset22","dataset23")){
       print(MutationType_pm(dataset_mu = datasets_mu[[ds]],gene = gene))
      }else if(ds %in% c("dataset11","dataset13","dataset24")){
        print(TMB_pm(datasets[[ds]],gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      }
      dev.off()
      shinyjs::enable("ref_mut_pm_down")
    }
  )
  
  output$SPLOT_pm = renderReactable({
    ds = ref_pm_para()$datasource_pm
    gene = ref_pm_para()$gene_pm
    Mut_type_ref_pm = ref_pm_para()$Mut_type_ref_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(! ds %in% c("dataset11","dataset13"),message = FALSE))
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_datatable = datasets_mu[[ds]][datasets_mu[[ds]]$Hugo_Symbol %in% pathway_list[[gene]] & datasets_mu[[ds]]$Variant_Classification %in% Mut_type_ref_pm,]
    for(i in colnames(tmp_datatable)){
      if(!i %in% c("Start_Position","End_Position","Position","ID")){
        tmp_datatable[,i] = as.factor(tmp_datatable[,i])
      }
    }

      tmp = reactable( tmp_datatable,
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
                                              align = "center",
                                              cell = function(value) format(value,digits = 2))
                       
      )
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)
  }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_cohort_pm(),width3(),ref_pm_para()$Mut_type_ref_pm)
  
  output$ref_mut_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      paste0(gene,"_",input$dataset_pm,".","mut", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("ref_mut_pm_tabdown")
      ds = ref_pm_para()$datasource_pm
      gene = ref_pm_para()$gene_pm
      Mut_type_ref_pm = ref_pm_para()$Mut_type_ref_pm
      

      write.csv(x = datasets_mu[[ds]][datasets_mu[[ds]]$Hugo_Symbol %in% pathway_list[[gene]] & datasets_mu[[ds]]$Variant_Classification %in% Mut_type_ref_pm,],file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("ref_mut_pm_tabdown")
    }
  )
  
  ###################################################ref DEG ################################################
  
  output$DEG_pm_uitabledown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    downloadButton('ref_diff_pm_tabdown',label = 'Download Table')
  })
  
  observe({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm

    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_id4_pm")
      createAlert(session, "warning4_pm", "warning_id4_pm", title = "Warning",style = "danger",
                  content = paste("The dataset",dataset_name2[[ds]],"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_id4_pm")
    }
    
    if(!any(pathway_list[[gene]] %in% colnames(datasets[[ds]]))){
      closeAlert(session,"warning_id4_pm")
      createAlert(session, "warning4_pm", "warning_id4_pm", title = "Warning",style = "danger",
                  content = paste("The genes of",gene,"does not exist in",dataset_name2[[ds]],sep = " "), append = FALSE)
    }

  })
  
  observe({
    
    gene = ref_pm_para()$gene_pm
    ref_cohort_pm = ref_cohort_pm()
    
    if(length(ref_cohort_pm$mut) <3 | length(ref_cohort_pm$wt) < 3){
      closeAlert(session,"warning_id4_pm")
      createAlert(session, "warning4_pm", "warning_id4_pm", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }
  })
  
  DEG_table = reactive({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    min.pct_ref_pm = min.pct_ref_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
   
    tmp = DEG_ref_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,min.pct = min.pct_ref_pm,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)
    
    }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_cohort_pm(),min.pct_ref_pm())
  
  output$DEG_tab_ref_pm = renderReactable({
    DEG_table = DEG_table()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table,3),
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
                                            )
                     
    )
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(DEG_table())
  
  output$ref_diff_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      paste0(gene,"_",input$dataset_pm,".","diff", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("ref_diff_pm_tabdown")
      write.csv(x = DEG_table(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("ref_diff_pm_tabdown")
    }
  )
  
  
  output$volcano_ref_pm = renderPlotly({
    ds = ref_pm_para()$datasource_pm
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    tab = DEG_table()
    FC_ref_pm = FC_ref_pm()
    pvalue_ref_pm = pvalue_ref_pm()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_ref_pm,pvalue = pvalue_ref_pm)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)
  }) %>% bindCache(ref_pm_para()$datasource_pm,DEG_table(),FC_ref_pm(),pvalue_ref_pm())
  
  #################################################ref GSEA######################################################
  output$GSEA_pm_uitabledown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    downloadButton('ref_gsea_pm_tabdown',label = 'Download Table')
  })
  
  observe({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    
    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_id5_pm")
      createAlert(session, "warning5_pm", "warning_id5_pm", title = "Warning",style = "danger",
                  content = paste("The dataset",dataset_name2[[ds]],"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_id5_pm")
    }
    
    if(!any(pathway_list[[gene]] %in% colnames(datasets[[ds]]))){
      closeAlert(session,"warning_id5_pm")
      createAlert(session, "warning5_pm", "warning_id5_pm", title = "Warning",style = "danger",
                  content = paste("The genes of",gene,"does not exist in",dataset_name2[[ds]],sep = " "), append = FALSE)
    }
    
  })
  
  observe({
    
    gene = ref_pm_para()$gene_pm
    ref_cohort_pm = ref_cohort_pm()
    
    if(length(ref_cohort_pm$mut) <3 | length(ref_cohort_pm$wt) < 3){
      closeAlert(session,"warning_id5_pm")
      createAlert(session, "warning5_pm", "warning_id5_pm", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }
  })
  
  observeEvent(input$ref_pm_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "ref_gsea_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "ref_gsea_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "ref_gsea_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "ref_gsea_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "ref_gsea_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_ref_pm"),type = 1)
    ))
  })
  
  observeEvent(input$ref_gsea_pm_file,{
    if(input$ref_gsea_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_gsea_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_gsea_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_gsea_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_gsea_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_gsea_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_gsea_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  status_file <- tempfile()
  write("", status_file)
  
  
  tmp_gsea_pm = reactive({
    diff = DEG_table()
    FC = diff$logFC
    names(FC) = rownames(diff)
    FC = sort(FC,decreasing = T)
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::disable("pathway_gsea_ref_single");shinyjs::disable("pathway_gsea_tcga");shinyjs::disable("pathway_gsea_cptac_rna");shinyjs::disable("pathway_gsea_cptac_protein");shinyjs::disable("pathway_gsea_ref_pm");shinyjs::disable("pathway_gsea_tcga_pm");shinyjs::disable("pathway_gsea_cptac_rna_pm");shinyjs::disable("pathway_gsea_cptac_protein_pm");shinyjs::disable("pathway_gsea_tcga_subtype");shinyjs::disable("pathway_gsea_cptac_rna_subtype");shinyjs::disable("pathway_gsea_cptac_protein_subtype");
    
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    if( length(GSEA_use) < 5 | useid %in% names(GSEA_use) ){
      
      installr::kill_pid(pid = scan(status_file))
      tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_ref_pm]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
      write(tmp$get_pid(), status_file)
      
    }else{
      
      closeAlert(session,"warning_id5_pm")
      createAlert(session, "warning5_pm", "warning_id5_pm", title = "Warning",style = "danger",
                  content = "Sorry, GSEA function is temporarily unavailable due to high user traffic. Please try again later.", append = FALSE)
      
      tmp = NULL
      
    }

    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::enable("pathway_gsea_ref_single");shinyjs::enable("pathway_gsea_tcga");shinyjs::enable("pathway_gsea_cptac_rna");shinyjs::enable("pathway_gsea_cptac_protein");shinyjs::enable("pathway_gsea_ref_pm");shinyjs::enable("pathway_gsea_tcga_pm");shinyjs::enable("pathway_gsea_cptac_rna_pm");shinyjs::enable("pathway_gsea_cptac_protein_pm");shinyjs::enable("pathway_gsea_tcga_subtype");shinyjs::enable("pathway_gsea_cptac_rna_subtype");shinyjs::enable("pathway_gsea_cptac_protein_subtype");
    
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  })
  
  one_loader = reactiveVal(0)
  output$GSEA_tab_ref_pm = renderReactable({
    
    if(is.null(tmp_gsea_pm())){return(NULL)}
    if(tmp_gsea_pm()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader() == 0){
        
        waiter_show(id = "ref_pm_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader(one_loader()+1)
        GSEA_use[[useid]] <<- TRUE
        return(NULL)
        
        
      }else if(one_loader() == 1){
        
        shinyjs::disable(id = "ref_gsea_pm_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "ref_pm_GSEA_loader")
      
      tmp_gsea_pm = tmp_gsea_pm()$get_result()@result
      tmp_gsea_pm$ID = as.factor(tmp_gsea_pm$ID)
      tmp_gsea_pm$Description = as.factor(tmp_gsea_pm$Description)
      tmp_gsea_pm$Plot = NA
      tmp_gsea_pm = tmp_gsea_pm[,c(ncol(tmp_gsea_pm),3:(ncol(tmp_gsea_pm)-2))]
      tmp = reactable( tmp_gsea_pm,
                       searchable = TRUE,
                       paginationType = "jump",
                       resizable = TRUE,
                       showPageSizeOptions = FALSE, 
                       highlight = TRUE,
                       bordered = TRUE,
                       outlined = TRUE,
                       striped = TRUE,
                       defaultColDef = colDef(minWidth = 150,headerStyle = list(background = "#033c73",color = "white"),
                                              sortNALast = TRUE,
                                              align = "center",
                                              cell = function(value) format(value,digits = 2)),
                       columns = list(
                         Plot=colDef(cell = function() htmltools::tags$button("GSEA Plot"))
                       ),
                       onClick = JS(
                         "function(rowInfo, column) {
                         if (column.id !== 'Plot') {return}
                         //window.alert('Details for row ' + rowInfo.index + ':\\n' + JSON.stringify(rowInfo.values, null, 2))
                         if (window.Shiny) {
                           Shiny.setInputValue('ref_pm_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      
      shinyjs::enable(id = "ref_gsea_pm_tabdown")
      one_loader(0)
      write("", status_file)
      
      GSEA_use[[useid]] <<- NULL
      return(tmp)
    }
    
  })
  
  output$ref_gsea_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      paste0(gene,"_",input$dataset_pm,".","gsea", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("ref_gsea_pm_tabdown")
      write.csv(x = tmp_gsea_pm()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("ref_gsea_pm_tabdown")
    }
  )
  
  
  row_num_pm = eventReactive(input$ref_pm_GSEA_show,{input$ref_pm_GSEA_show$index})%>% debounce(500)
  
  output$GSEA_plot_ref_pm = renderPlot({
    my_gseaplot2(tmp_gsea_pm()$get_result(),geneSetID = tmp_gsea_pm()$get_result()@result$ID[row_num_pm()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_pm()$get_result()@result$ID[row_num_pm()])
  })
  
  output$ref_gsea_pm_down = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      
      if(input$ref_gsea_pm_file == 'pdf'){
        paste0(gene,"_",input$dataset_pm,".","gsea", '.',"pdf")
      }else if(input$ref_gsea_pm_file == 'png'){
        paste0(gene,"_",input$dataset_pm,".","gsea", '.',"png")
      }else if(input$ref_gsea_pm_file == 'jpeg'){
        paste0(gene,"_",input$dataset_pm,".","gsea", '.',"jpeg")
      }else{
        paste0(gene,"_",input$dataset_pm,".","gsea", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("ref_gsea_pm_down")
      gene = ref_pm_para()$gene_pm
      ds = ref_pm_para()$datasource_pm
      ref_cohort_pm = ref_cohort_pm()
      
      if(input$ref_gsea_pm_res <= 300){
        r = input$ref_gsea_pm_res
      }else{
        r = 300
      }
      if(input$ref_gsea_pm_file == "pdf" & input$ref_gsea_pm_height <= 50){
        h = input$ref_gsea_pm_height
      }else if(input$ref_gsea_pm_file == "pdf" & input$ref_gsea_pm_height > 50){
        h = 50
      }else if(input$ref_gsea_pm_file != "pdf" & input$ref_gsea_pm_height <= 3000){
        h = input$ref_gsea_pm_height
      }else{
        h = 3000
      }
      if(input$ref_gsea_pm_file == "pdf" & input$ref_gsea_pm_width <= 50){
        w = input$ref_gsea_pm_width
      }else if(input$ref_gsea_pm_file == "pdf" & input$ref_gsea_pm_width > 50){
        w = 50
      }else if(input$ref_gsea_pm_file != "pdf" & input$ref_gsea_pm_width <= 3000){
        w = input$ref_gsea_pm_width
      }else{
        w = 3000
      }
      
      if(input$ref_gsea_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_gsea_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_gsea_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(
        my_gseaplot2(tmp_gsea_pm()$get_result(),geneSetID = tmp_gsea_pm()$get_result()@result$ID[row_num_pm()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_pm()$get_result()@result$ID[row_num_pm()])
      )
      
      
      dev.off()
      shinyjs::enable("ref_gsea_pm_down")
    }
  )
  ################################################ref immune #################################################
  output$inf_pm_uidown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_inf_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_inf_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_inf_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_inf_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_inf_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$ref_inf_pm_file,{
    if(input$ref_inf_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_inf_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_inf_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_inf_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_inf_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_inf_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_inf_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$sig_pm_uidown = renderUI({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3 ,message = FALSE))
    validate(need(any(pathway_list[[gene]] %in% colnames(datasets[[ds]])),message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_sig_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_sig_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sig_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sig_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_sig_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$ref_sig_pm_file,{
    if(input$ref_sig_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_sig_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_sig_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_sig_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_sig_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_sig_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_sig_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    
    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_id6_pm")
      createAlert(session, "warning6_pm", "warning_id6_pm", title = "Warning",style = "danger",
                  content = paste("The dataset",dataset_name2[[ds]],"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_id6_pm")
    }
    
    if(!any(pathway_list[[gene]] %in% colnames(datasets[[ds]]))){
      closeAlert(session,"warning_id6_pm")
      createAlert(session, "warning6_pm", "warning_id6_pm", title = "Warning",style = "danger",
                  content = paste("The genes of",gene,"does not exist in",dataset_name2[[ds]],sep = " "), append = FALSE)
    }
    
  })
  
  observe({
    
    gene = ref_pm_para()$gene_pm
    ref_cohort_pm = ref_cohort_pm()
    
    if(length(ref_cohort_pm$mut) <3 | length(ref_cohort_pm$wt) < 3){
      closeAlert(session,"warning_id6_pm")
      createAlert(session, "warning6_pm", "warning_id6_pm", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }
  })
  
  output$immune_infiltration_datasets_all_pm = renderPlotly({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_all_datasets_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_infiltration",dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    tmp = ggplotly(p) %>% 
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0), height = 800) %>% 
            rangeslider(start = 0,end = 15.5) %>% plotly::config(displayModeBar = FALSE)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    
    return(tmp)
  }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_cohort_pm())
  
  output$ref_inf_pm_down = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      
      if(input$ref_inf_pm_file == 'pdf'){
        paste0(gene,"_",input$dataset_pm,".","inf", '.',"pdf")
      }else if(input$ref_inf_pm_file == 'png'){
        paste0(gene,"_",input$dataset_pm,".","inf", '.',"png")
      }else if(input$ref_inf_pm_file == 'jpeg'){
        paste0(gene,"_",input$dataset_pm,".","inf", '.',"jpeg")
      }else{
        paste0(gene,"_",input$dataset_pm,".","inf", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("ref_inf_pm_down")
      gene = ref_pm_para()$gene_pm
      ds = ref_pm_para()$datasource_pm
      ref_cohort_pm = ref_cohort_pm()
      
      if(input$ref_inf_pm_res <= 300){
        r = input$ref_inf_pm_res
      }else{
        r = 300
      }
      if(input$ref_inf_pm_file == "pdf" & input$ref_inf_pm_height <= 50){
        h = input$ref_inf_pm_height
      }else if(input$ref_inf_pm_file == "pdf" & input$ref_inf_pm_height > 50){
        h = 50
      }else if(input$ref_inf_pm_file != "pdf" & input$ref_inf_pm_height <= 3000){
        h = input$ref_inf_pm_height
      }else{
        h = 3000
      }
      if(input$ref_inf_pm_file == "pdf" & input$ref_inf_pm_width <= 50){
        w = input$ref_inf_pm_width
      }else if(input$ref_inf_pm_file == "pdf" & input$ref_inf_pm_width > 50){
        w = 50
      }else if(input$ref_inf_pm_file != "pdf" & input$ref_inf_pm_width <= 3000){
        w = input$ref_inf_pm_width
      }else{
        w = 3000
      }
      
      if(input$ref_inf_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_inf_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_inf_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(immune_all_datasets_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_infiltration",dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      
      dev.off()
      shinyjs::enable("ref_inf_pm_down")
    }
  )
  
  output$immune_signature_datasets_all_pm = renderPlotly({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_all_datasets_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_pathway",dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    tmp = ggplotly(p) %>% 
      layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0), height = 800) %>% 
      rangeslider(start = 0,end = nrow(datasets_rna_wes[[ds]]$immune_pathway)+1) %>% plotly::config(displayModeBar = FALSE)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_cohort_pm())
  
  
  output$ref_sig_pm_down = downloadHandler(
    
    filename = function(){
      gene = ref_pm_para()$gene_pm
      
      if(input$ref_sig_pm_file == 'pdf'){
        paste0(gene,"_",input$dataset_pm,".","sig", '.',"pdf")
      }else if(input$ref_sig_pm_file == 'png'){
        paste0(gene,"_",input$dataset_pm,".","sig", '.',"png")
      }else if(input$ref_sig_pm_file == 'jpeg'){
        paste0(gene,"_",input$dataset_pm,".","sig", '.',"jpeg")
      }else{
        paste0(gene,"_",input$dataset_pm,".","sig", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("ref_sig_pm_down")
      gene = ref_pm_para()$gene_pm
      ds = ref_pm_para()$datasource_pm
      ref_cohort_pm = ref_cohort_pm()
      
      if(input$ref_sig_pm_res <= 300){
        r = input$ref_sig_pm_res
      }else{
        r = 300
      }
      if(input$ref_sig_pm_file == "pdf" & input$ref_sig_pm_height <= 50){
        h = input$ref_sig_pm_height
      }else if(input$ref_sig_pm_file == "pdf" & input$ref_sig_pm_height > 50){
        h = 50
      }else if(input$ref_sig_pm_file != "pdf" & input$ref_sig_pm_height <= 3000){
        h = input$ref_sig_pm_height
      }else{
        h = 3000
      }
      if(input$ref_sig_pm_file == "pdf" & input$ref_sig_pm_width <= 50){
        w = input$ref_sig_pm_width
      }else if(input$ref_sig_pm_file == "pdf" & input$ref_sig_pm_width > 50){
        w = 50
      }else if(input$ref_sig_pm_file != "pdf" & input$ref_sig_pm_width <= 3000){
        w = input$ref_sig_pm_width
      }else{
        w = 3000
      }
      
      if(input$ref_sig_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_sig_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_sig_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(immune_all_datasets_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_pathway",dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt))
      
      dev.off()
      shinyjs::enable("ref_sig_pm_down")
    }
  )
  
  output$immune_infiltration_datasets_one_pm = renderPlot({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_datasets_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_infiltration",selection = input$immune_infiltration_datasets_input_pm,outlier = input$outlier_datasets_pm,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_cohort_pm(),input$immune_infiltration_datasets_input_pm,input$outlier_datasets_pm)
  
  output$immune_signature_datasets_one_pm = renderPlot({
    gene = ref_pm_para()$gene_pm
    ds = ref_pm_para()$datasource_pm
    ref_cohort_pm = ref_cohort_pm()
    validate(need(length(ref_cohort_pm$mut) >=3 & length(ref_cohort_pm$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_datasets_pm(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_pathway",selection = input$immune_signature_datasets_input_pm,outlier = input$outlier_datasets_pm,dataset_mu = datasets_mu[[ds]],mut = ref_cohort_pm$mut,wt = ref_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(ref_pm_para()$datasource_pm,ref_pm_para()$gene_pm,ref_cohort_pm(),input$immune_signature_datasets_input_pm,input$outlier_datasets_pm)
}