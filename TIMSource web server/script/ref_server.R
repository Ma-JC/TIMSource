Ref_datasets_server <- function(input,output,session,ref_single_para,datasets,datasets_mu,pathway_database,min.pct_ref_single,FC_ref_single,pvalue_ref_single,useid){
  
  useid <- paste(useid,"RefSingle",sep = "_")
  
  width = reactive({if(ref_single_para()$datasource %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset4","dataset6","dataset9","dataset12","dataset14","dataset15","dataset17","dataset18","dataset19","dataset20","dataset21","dataset22","dataset23","dataset24")){900}else{500}})
  width2 = reactive({if(ref_single_para()$datasource %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset7","dataset12","dataset14","dataset16","dataset24")){1200}else{600}})
  width3 = reactive({if(ref_single_para()$datasource %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9","dataset3.1","dataset3.2","dataset3.3","dataset3",
                                            "dataset1","dataset2","dataset4","dataset5","dataset6","dataset7","dataset8","dataset9","dataset10","dataset12","dataset14","dataset15","dataset16","dataset17","dataset20","dataset21")){1600}else{1000}})
  
  ref_cohort = reactive({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    Mut_type_ref_single = ref_single_para()$Mut_type_ref_single
    Wild_type_ref_single = ref_single_para()$Wild_type_ref_single
    Therapy_type_ref_single = ref_single_para()$Therapy_type_ref_single
   
    ref_cohort_cal(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],Mut_type = Mut_type_ref_single,Wild_type = Wild_type_ref_single,therapy_type=Therapy_type_ref_single)
    }) %>% bindCache(ref_single_para()$datasource,ref_single_para()$gene,ref_single_para()$Mut_type_ref_single,ref_single_para()$Wild_type_ref_single,ref_single_para()$Therapy_type_ref_single)
  
  ############################Survival#########################################
  output$sur_single_uidown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(!ds %in% c(),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_sur_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_sur_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sur_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sur_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_sur_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$ref_sur_single_file,{
    if(input$ref_sur_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_sur_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_sur_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_sur_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_sur_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_sur_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_sur_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    if(ds %in% c("dataset5","dataset7","dataset13","dataset16")){
      closeAlert(session,"warning_id")
      createAlert(session, "warning", "warning_id", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide overall survival"), append = FALSE)
    }else if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                       "dataset1","dataset8","dataset10","dataset11")){
      closeAlert(session,"warning_id")
      createAlert(session, "warning", "warning_id", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide progression-free survival"), append = FALSE)
    }else{
      closeAlert(session,"warning_id")
    }

    if(length(ref_cohort$mut) <3 | length(ref_cohort$wt) < 3){
      closeAlert(session,"warning_id")
      createAlert(session, "warning", "warning_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }
  })
  
  output$Survival = renderPlot({
    
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
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
      p1 = os_survival(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
      p2 = pfs_survival(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
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
      tmp = os_survival(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
    }else if(ds %in%  c("dataset5","dataset7","dataset13","dataset16")){
      tmp = pfs_survival(dataset = datasets[[ds]],gene = gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
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
  

  output$ref_sur_single_down = downloadHandler(

    filename = function(){
      if(input$ref_sur_single_file == 'pdf'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","sur", '.',"pdf")
      }else if(input$ref_sur_single_file == 'png'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","sur", '.',"png")
      }else if(input$ref_sur_single_file == 'jpeg'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","sur", '.',"jpeg")
      }else{
        paste0(ref_single_para()$gene,"_",input$dataset,".","sur", '.',"tiff")
      }

    },
    content = function(file){

      shinyjs::disable("ref_sur_single_down")
      if(input$ref_sur_single_res <= 300){
        r = input$ref_sur_single_res
      }else{
        r = 300
      }
      if(input$ref_sur_single_file == "pdf" & input$ref_sur_single_height <= 50){
        h = input$ref_sur_single_height
      }else if(input$ref_sur_single_file == "pdf" & input$ref_sur_single_height > 50){
        h = 50
      }else if(input$ref_sur_single_file != "pdf" & input$ref_sur_single_height <= 3000){
        h = input$ref_sur_single_height
      }else{
        h = 3000
      }
      if(input$ref_sur_single_file == "pdf" & input$ref_sur_single_width <= 50){
        w = input$ref_sur_single_width
      }else if(input$ref_sur_single_file == "pdf" & input$ref_sur_single_width > 50){
        w = 50
      }else if(input$ref_sur_single_file != "pdf" & input$ref_sur_single_width <= 3000){
        w = input$ref_sur_single_width
      }else{
        w = 3000
      }

      if(input$ref_sur_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_sur_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_sur_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      if(ref_single_para()$datasource %in%  c("dataset2","dataset3","dataset4","dataset6","dataset9","dataset12","dataset14","dataset15","dataset17","dataset18","dataset19","dataset20","dataset21","dataset22","dataset23","dataset24")){
        p1 = os_survival(dataset = datasets[[ref_single_para()$datasource]],gene = ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
        p2 = pfs_survival(dataset = datasets[[ref_single_para()$datasource]],gene = ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
        layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
        print((p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout))
      }else if(ref_single_para()$datasource %in%  c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                                    "dataset1","dataset8","dataset10","dataset11")){
        print(os_survival(dataset = datasets[[ref_single_para()$datasource]],gene = ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt))
      }else if(ref_single_para()$datasource %in%  c("dataset5","dataset7","dataset13","dataset16")){
        print(pfs_survival(dataset = datasets[[ref_single_para()$datasource]],gene = ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt))
      }
      dev.off()

      shinyjs::enable("ref_sur_single_down")

    }
  )

  # output$warning = renderText({
  #   gene = ref_single_para()$gene
  #   ds = ref_single_para()$datasource
  #   if(ds %in% c("dataset5","dataset7","dataset13","dataset16")){
  #     return(paste(ds,"didn't provide overall survival"))
  #   }else if(ds %in% c("dataset1","dataset8","dataset10","dataset11")){
  #     return(paste(ds,"didn't provide progression-free survival"))
  #   }
  #
  # })


  ##########################################Drug response####################################
  output$res_single_uidown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(!ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                             "dataset1","dataset13"),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_res_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_res_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_res_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_res_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_res_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })

  observeEvent(input$ref_res_single_file,{
    if(input$ref_res_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_res_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_res_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_res_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_res_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_res_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_res_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })

  observe({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    if(ds %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9",
                 "dataset1","dataset13")){
      closeAlert(session,"warning_id2")
      createAlert(session, "warning2", "warning_id2", title = "Warning",style = "danger",
                  content = paste(dataset_name2[[ds]],"did not provide any data on the efficacy of immunotherapy"), append = FALSE)

    }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
      closeAlert(session,"warning_id2")
      createAlert(session, "warning2", "warning_id2", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide data about RECIST"), append = FALSE)

    }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
      closeAlert(session,"warning_id2")
      createAlert(session, "warning2", "warning_id2", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide data about DCB"), append = FALSE)
    }else{
      closeAlert(session,"warning_id2")
    }

    if(length(ref_cohort$mut) <3 | length(ref_cohort$wt) < 3){
      closeAlert(session,"warning_id2")
      createAlert(session, "warning2", "warning_id2", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }

  })


  output$response = renderPlot({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    width2 = width2()
    ref_cohort = ref_cohort()
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
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
      p1 = CRPR(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
      p2 = DCB(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
      tmp = p1 + p2
    }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
      tmp = DCB(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
    }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
      tmp = CRPR(datasets[[ds]],gene,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
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

  output$ref_res_single_down = downloadHandler(
    filename = function(){
      if(input$ref_res_single_file == 'pdf'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","res", '.',"pdf")
      }else if(input$ref_res_single_file == 'png'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","res", '.',"png")
      }else if(input$ref_res_single_file == 'jpeg'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","res", '.',"jpeg")
      }else{
        paste0(ref_single_para()$gene,"_",input$dataset,".","res", '.',"tiff")
      }

    },
    content = function(file){
      shinyjs::disable("ref_res_single_down")
      if(input$ref_res_single_res <= 300){
        r = input$ref_res_single_res
      }else{
        r = 300
      }
      if(input$ref_res_single_file == "pdf" & input$ref_res_single_height <= 50){
        h = input$ref_res_single_height
      }else if(input$ref_res_single_file == "pdf" & input$ref_res_single_height > 50){
        h = 50
      }else if(input$ref_res_single_file != "pdf" & input$ref_res_single_height <= 3000){
        h = input$ref_res_single_height
      }else{
        h = 3000
      }
      if(input$ref_res_single_file == "pdf" & input$ref_res_single_width <= 50){
        w = input$ref_res_single_width
      }else if(input$ref_res_single_file == "pdf" & input$ref_res_single_width > 50){
        w = 50
      }else if(input$ref_res_single_file != "pdf" & input$ref_res_single_width <= 3000){
        w = input$ref_res_single_width
      }else{
        w = 3000
      }

      if(input$ref_res_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_res_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_res_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }


      if(ref_single_para()$datasource %in% c("dataset2","dataset3.1","dataset3.2","dataset3.3","dataset3","dataset7","dataset12","dataset14","dataset16","dataset24")){
        p1 = CRPR(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
        p2 = DCB(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
        print(p1 + p2)
      }else if(ref_single_para()$datasource %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
        print(DCB(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt))
      }else if(ref_single_para()$datasource %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20","dataset21","dataset22","dataset23")){
        print(CRPR(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt))
      }

      dev.off()
      shinyjs::enable("ref_res_single_down")
    }
  )
  #
  # output$warning2 = renderText({
  #   gene = ref_single_para()$gene
  #   ds = ref_single_para()$datasource
  #   if(ds %in% c("dataset1","dataset13")){
  #     paste(ds,"did not provide any data on the efficacy of immunotherapy")
  #   }else if(ds %in% c("dataset5","dataset6","dataset9","dataset10","dataset15","dataset19")){
  #     paste(ds,"didn't provide data about RECIST")
  #   }else if(ds %in% c("dataset4","dataset8","dataset11","dataset17","dataset18","dataset20")){
  #     paste(ds,"didn't provide data about DCB")
  #   }
  #
  # })

  #############################################Mutation##############################################
  output$mut_single_uidown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(!ds %in% c(),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_mut_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_mut_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_mut_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_mut_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_mut_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })

  observeEvent(input$ref_mut_single_file,{
    if(input$ref_mut_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_mut_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_mut_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_mut_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_mut_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_mut_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_mut_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })

  output$mut_single_uitabledown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(!ds %in% c("dataset11","dataset13"),message = FALSE))
    downloadButton('ref_mut_single_tabdown',label = 'Download Table')
  })

  observe({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    if(ds %in% c("dataset18","dataset19","dataset22","dataset23")){
      closeAlert(session,"warning_id3")
      createAlert(session, "warning3", "warning_id3", title = "Tips",
                  content = paste(dataset_name2[[ds]],"did not provide TMB data"), append = FALSE)

    }else if(ds %in% c("dataset11","dataset13")){
      closeAlert(session,"warning_id3")
      createAlert(session, "warning3", "warning_id3", title = "Tips",
                  content = paste(dataset_name2[[ds]],"didn't provide detial information about gene mutation"), append = FALSE)

    }else{
      closeAlert(session,"warning_id3")
    }

    if(length(ref_cohort$mut) <3 | length(ref_cohort$wt) < 3){
      closeAlert(session,"warning_id3")
      createAlert(session, "warning3", "warning_id3", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }

  })

  output$TMB = renderPlot({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    width3 = width3()
    ref_cohort = ref_cohort()
    outline = input$outline
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))

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
        p1 = TMB(datasets[[ds]],gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
      }else{p1 = TMB(datasets[[ds]],gene,remove_outlier = FALSE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)}

      p2 = MutationType(dataset_mu = datasets_mu[[ds]],gene = gene)
      tmp = p1 + p2
    }else if(ds %in% c("dataset18","dataset19","dataset22","dataset23")){
      tmp = MutationType(dataset_mu = datasets_mu[[ds]],gene = gene)
    }else if(ds %in% c("dataset11","dataset13","dataset24")){
      tmp = TMB(datasets[[ds]],gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
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



  output$ref_mut_single_down = downloadHandler(
    filename = function(){
      if(input$ref_mut_single_file == 'pdf'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","mut", '.',"pdf")
      }else if(input$ref_mut_single_file == 'png'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","mut", '.',"png")
      }else if(input$ref_mut_single_file == 'jpeg'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","mut", '.',"jpeg")
      }else{
        paste0(ref_single_para()$gene,"_",input$dataset,".","mut", '.',"tiff")
      }

    },
    content = function(file){
      shinyjs::disable("ref_mut_single_down")
      if(input$ref_mut_single_res <= 300){
        r = input$ref_mut_single_res
      }else{
        r = 300
      }
      if(input$ref_mut_single_file == "pdf" & input$ref_mut_single_height <= 50){
        h = input$ref_mut_single_height
      }else if(input$ref_mut_single_file == "pdf" & input$ref_mut_single_height > 50){
        h = 50
      }else if(input$ref_mut_single_file != "pdf" & input$ref_mut_single_height <= 3000){
        h = input$ref_mut_single_height
      }else{
        h = 3000
      }
      if(input$ref_mut_single_file == "pdf" & input$ref_mut_single_width <= 50){
        w = input$ref_mut_single_width
      }else if(input$ref_mut_single_file == "pdf" & input$ref_mut_single_width > 50){
        w = 50
      }else if(input$ref_mut_single_file != "pdf" & input$ref_mut_single_width <= 3000){
        w = input$ref_mut_single_width
      }else{
        w = 3000
      }

      if(input$ref_mut_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_mut_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_mut_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      if(ref_single_para()$datasource %in% c("dataset1.1","dataset1.2","dataset1.3","dataset1.4","dataset1.5","dataset1.6","dataset1.7","dataset1.8","dataset1.9","dataset3.1","dataset3.2","dataset3.3","dataset3",
                             "dataset1","dataset2","dataset4","dataset5","dataset6","dataset7","dataset8","dataset9","dataset10","dataset12","dataset14","dataset15","dataset16","dataset17","dataset20","dataset21")){

        if(input$outline){
          p1 = TMB(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
        }else{p1 = TMB(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,remove_outlier = FALSE,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)}

        p2 = MutationType(dataset_mu = datasets_mu[[ref_single_para()$datasource]],gene = ref_single_para()$gene)
        print(p1 + p2)
      }else if(ref_single_para()$datasource %in% c("dataset18","dataset19","dataset22","dataset23")){
        print(MutationType(dataset_mu = datasets_mu[[ref_single_para()$datasource]],gene = ref_single_para()$gene))
      }else if(ref_single_para()$datasource %in% c("dataset11","dataset13","dataset24")){
        print(TMB(datasets[[ref_single_para()$datasource]],ref_single_para()$gene,remove_outlier = TRUE,dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt))
      }

      dev.off()
      shinyjs::enable("ref_mut_single_down")
    }
  )


  output$SPLOT = renderReactable({
    ds = ref_single_para()$datasource
    gene = ref_single_para()$gene
    Mut_type_ref_single = ref_single_para()$Mut_type_ref_single
    ref_cohort = ref_cohort()
    validate(need(! ds %in% c("dataset11","dataset13"),message = FALSE))
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))

    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");

      tmp_datatable = datasets_mu[[ds]][datasets_mu[[ds]]$Hugo_Symbol %in% gene & datasets_mu[[ds]]$Variant_Classification %in% Mut_type_ref_single,]
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

  }) %>% bindCache(ref_single_para()$datasource,ref_single_para()$gene,ref_cohort(),ref_single_para()$Mut_type_ref_single)


  output$ref_mut_single_tabdown = downloadHandler(

    filename = function(){
      gene = ref_single_para()$gene
      paste0(gene,"_",input$dataset,".","mut", '.',"csv")

    },
    content = function(file){
      shinyjs::disable("ref_mut_single_tabdown")
      ds = ref_single_para()$datasource
      gene = ref_single_para()$gene
      Mut_type_ref_single = ref_single_para()$Mut_type_ref_single
      ref_cohort = ref_cohort()


      write.csv(x = datasets_mu[[ds]][datasets_mu[[ds]]$Hugo_Symbol %in% gene & datasets_mu[[ds]]$Variant_Classification %in% Mut_type_ref_single,],file = file, sep = ',', col.names = T, row.names = T, quote = F)

      shinyjs::enable("ref_mut_single_tabdown")
      }
  )
  # observeEvent(input$SPLOT_rows_selected,{message(input$SPLOT_rows_selected)})

  ###################################################ref DEG ################################################

  observe({
    ds = ref_single_para()$datasource
    gene = ref_single_para()$gene
    ref_cohort = ref_cohort()
    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_id4")
      createAlert(session, "warning4", "warning_id4", title = "Warning",style = "danger",
                  content = paste(dataset_name2[[ds]],"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_id4")
    }

    if(length(ref_cohort$mut) <3 | length(ref_cohort$wt) < 3){
      closeAlert(session,"warning_id4")
      createAlert(session, "warning4", "warning_id4", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }

  })

  output$DEG_single_uitabledown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    downloadButton('ref_diff_single_tabdown',label = 'Download Table')
  })

  DEG_table = reactive({
    ds = ref_single_para()$datasource
    gene = ref_single_para()$gene
    ref_cohort = ref_cohort()
    min.pct_ref_single = min.pct_ref_single()
    # validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " ")))
    # validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),paste("The dataset",ds,"didn't provide RNA-seq data",sep = " ")))
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");

    tmp = DEG_ref(datasets_rna_wes = datasets_rna_wes,dataset = ds,min.pct = min.pct_ref_single,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)

    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");

    return(tmp)


    }) %>% bindCache(ref_single_para()$datasource,ref_single_para()$gene,ref_cohort(),min.pct_ref_single())

  output$DEG_tab_ref_single = renderReactable({
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

  output$ref_diff_single_tabdown = downloadHandler(

    filename = function(){
      gene = ref_single_para()$gene
      paste0(gene,"_",input$dataset,".","diff", '.',"csv")

    },
    content = function(file){
      shinyjs::disable("ref_diff_single_tabdown")
      write.csv(x = DEG_table(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("ref_diff_single_tabdown")
    }
  )


  output$volcano_ref_single = renderPlotly({
    tab = DEG_table()
    FC_ref_single = FC_ref_single()
    pvalue_ref_single = pvalue_ref_single()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");

    tmp = volcano_plot(tab = tab,FC = FC_ref_single,pvalue = pvalue_ref_single)

    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");

    return(tmp)
  }) %>% bindCache(DEG_table(),FC_ref_single(),pvalue_ref_single())



  #################################################ref GSEA######################################################
  observe({
    ds = ref_single_para()$datasource
    gene = ref_single_para()$gene
    ref_cohort = ref_cohort()
    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_id5")
      createAlert(session, "warning5", "warning_id5", title = "Warning",style = "danger",
                  content = paste(dataset_name2[[ds]],"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_id5")
    }

    if(length(ref_cohort$mut) <3 | length(ref_cohort$wt) < 3){
      closeAlert(session,"warning_id5")
      createAlert(session, "warning5", "warning_id5", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }

  })

  output$GSEA_single_uitabledown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource

    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    downloadButton('ref_gsea_single_tabdown',label = 'Download Table')

  })


  status_file <- tempfile()
  write("", status_file)

  tmp_gsea = reactive({

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
      tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_ref_single]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
      write(tmp$get_pid(), status_file)
    }else{

      closeAlert(session,"warning_id5")
      createAlert(session, "warning5", "warning_id5", title = "Warning",style = "danger",
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
  output$GSEA_tab_ref_single = renderReactable({

    if(is.null(tmp_gsea())){return(NULL)}
    
    if(tmp_gsea()$is_alive()){

      invalidateLater(millis = 1000, session = session)

      if(one_loader() == 0){

        waiter_show(id = "ref_single_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader(one_loader()+1)
        
        GSEA_use[[useid]] <<- TRUE
        
        
        return(NULL)


        }else if(one_loader() == 1){

          shinyjs::disable(id = "ref_gsea_single_tabdown")

          }

      }else{
      waiter_hide(id = "ref_single_GSEA_loader")

      tmp_gsea = tmp_gsea()$get_result()@result
      tmp_gsea$ID = as.factor(tmp_gsea$ID)
      tmp_gsea$Description = as.factor(tmp_gsea$Description)
      tmp_gsea$Plot = NA
      tmp_gsea = tmp_gsea[,c(ncol(tmp_gsea),3:(ncol(tmp_gsea)-2))]
      tmp = reactable( tmp_gsea,defaultPageSize = 10,
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
                           Shiny.setInputValue('ref_single_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )

      )


      shinyjs::enable(id = "ref_gsea_single_tabdown")
      one_loader(0)
      write("", status_file)
      
      GSEA_use[[useid]] <<- NULL
      return(tmp)
    }


  })




  output$ref_gsea_single_tabdown = downloadHandler(

    filename = function(){
      gene = ref_single_para()$gene
      paste0(gene,"_",input$dataset,".","gsea", '.',"csv")

    },
    content = function(file){
      shinyjs::disable("ref_gsea_single_tabdown")
      write.csv(x = tmp_gsea()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("ref_gsea_single_tabdown")
    }
  )


  observeEvent(input$ref_single_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
               div(selectInput(inputId = "ref_gsea_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;"),
               div(numericInput(inputId = "ref_gsea_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
               div(numericInput(inputId = "ref_gsea_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
               div(numericInput(inputId = "ref_gsea_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
               div(downloadButton(outputId = "ref_gsea_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")

      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_ref_single"),type = 1)
    ))

  })

  observeEvent(input$ref_gsea_single_file,{
    if(input$ref_gsea_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_gsea_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_gsea_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_gsea_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_gsea_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_gsea_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_gsea_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })


  row_num = eventReactive(input$ref_single_GSEA_show,{input$ref_single_GSEA_show$index})%>% debounce(500)
  output$GSEA_plot_ref_single = renderPlot({
    my_gseaplot2(tmp_gsea()$get_result(),geneSetID = tmp_gsea()$get_result()@result$ID[row_num()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea()$get_result()@result$ID[row_num()])
  })


  output$ref_gsea_single_down = downloadHandler(
    filename = function(){
      if(input$ref_gsea_single_file == 'pdf'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","gsea", '.',"pdf")
      }else if(input$ref_gsea_single_file == 'png'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","gsea", '.',"png")
      }else if(input$ref_gsea_single_file == 'jpeg'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","gsea", '.',"jpeg")
      }else{
        paste0(ref_single_para()$gene,"_",input$dataset,".","gsea", '.',"tiff")
      }

    },
    content = function(file){
      shinyjs::disable("ref_gsea_single_down")
      if(input$ref_gsea_single_res <= 300){
        r = input$ref_gsea_single_res
      }else{
        r = 300
      }
      if(input$ref_gsea_single_file == "pdf" & input$ref_gsea_single_height <= 50){
        h = input$ref_gsea_single_height
      }else if(input$ref_gsea_single_file == "pdf" & input$ref_gsea_single_height > 50){
        h = 50
      }else if(input$ref_gsea_single_file != "pdf" & input$ref_gsea_single_height <= 3000){
        h = input$ref_gsea_single_height
      }else{
        h = 3000
      }
      if(input$ref_gsea_single_file == "pdf" & input$ref_gsea_single_width <= 50){
        w = input$ref_gsea_single_width
      }else if(input$ref_gsea_single_file == "pdf" & input$ref_gsea_single_width > 50){
        w = 50
      }else if(input$ref_gsea_single_file != "pdf" & input$ref_gsea_single_width <= 3000){
        w = input$ref_gsea_single_width
      }else{
        w = 3000
      }

      if(input$ref_gsea_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_gsea_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_gsea_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }


      print(    my_gseaplot2(tmp_gsea()$get_result(),geneSetID = tmp_gsea()$get_result()@result$ID[row_num()],
                             self.Description = "",
                             color="firebrick",
                             pvalue_table = T,
                             rel_heights = c(1, .2, 0.4),
                             title = tmp_gsea()$get_result()@result$ID[row_num()])
                )
      dev.off()
      shinyjs::enable("ref_gsea_single_down")
    }
  )

  ############################################## ref immune ##################################################
  observe({
    ds = ref_single_para()$datasource
    gene = ref_single_para()$gene
    ref_cohort = ref_cohort()
    if(!ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20")){
      closeAlert(session,"warning_id6")
      createAlert(session, "warning6", "warning_id6", title = "Warning",style = "danger",
                  content = paste(dataset_name2[[ds]],"didn't provide RNA-seq data",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_id6")
    }

    if(length(ref_cohort$mut) <3 | length(ref_cohort$wt) < 3){
      closeAlert(session,"warning_id6")
      createAlert(session, "warning6", "warning_id6", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype <3"), append = FALSE)
    }

  })

  output$inf_single_uidown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_inf_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_inf_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_inf_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_inf_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_inf_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })

  observeEvent(input$ref_inf_single_file,{
    if(input$ref_inf_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_inf_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_inf_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_inf_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_inf_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_inf_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_inf_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })

  output$sig_single_uidown = renderUI({
    ref_cohort = ref_cohort()
    ds = ref_single_para()$datasource
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))
    tagList(
      div(selectInput(inputId = "ref_sig_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "ref_sig_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sig_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "ref_sig_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "ref_sig_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })

  observeEvent(input$ref_sig_single_file,{
    if(input$ref_sig_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "ref_sig_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "ref_sig_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "ref_sig_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "ref_sig_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "ref_sig_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "ref_sig_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })

  output$immune_infiltration_datasets_all = renderPlotly({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    shinyjs::disable("sec")
    shinyjs::removeClass(id = "button_id",class = "fa fa-arrow-alt-circle-up")
    shinyjs::addClass(id = "button_id",class = "fa fa-spinner fa-spin")


    p = immune_all_datasets(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_infiltration",dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
    tmp = ggplotly(p) %>%
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0), height = 800) %>%
            rangeslider(start = 0,end = 15.5) %>% plotly::config(displayModeBar = FALSE)

    shinyjs::enable("sec")
    shinyjs::removeClass(id = "button_id",class = "fa fa-spinner fa-spin")
    shinyjs::addClass(id = "button_id",class = "fa fa-arrow-alt-circle-up")

    return(tmp)

  }) %>% bindCache(ref_single_para()$gene,ref_single_para()$datasource,ref_cohort())

  output$ref_inf_single_down = downloadHandler(
    filename = function(){
      if(input$ref_inf_single_file == 'pdf'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","inf", '.',"pdf")
      }else if(input$ref_inf_single_file == 'png'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","inf", '.',"png")
      }else if(input$ref_inf_single_file == 'jpeg'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","inf", '.',"jpeg")
      }else{
        paste0(ref_single_para()$gene,"_",input$dataset,".","inf", '.',"tiff")
      }

    },
    content = function(file){
      shinyjs::disable("ref_inf_single_down")
      if(input$ref_inf_single_res <= 300){
        r = input$ref_inf_single_res
      }else{
        r = 300
      }
      if(input$ref_inf_single_file == "pdf" & input$ref_inf_single_height <= 50){
        h = input$ref_inf_single_height
      }else if(input$ref_inf_single_file == "pdf" & input$ref_inf_single_height > 50){
        h = 50
      }else if(input$ref_inf_single_file != "pdf" & input$ref_inf_single_height <= 3000){
        h = input$ref_inf_single_height
      }else{
        h = 3000
      }
      if(input$ref_inf_single_file == "pdf" & input$ref_inf_single_width <= 50){
        w = input$ref_inf_single_width
      }else if(input$ref_inf_single_file == "pdf" & input$ref_inf_single_width > 50){
        w = 50
      }else if(input$ref_inf_single_file != "pdf" & input$ref_inf_single_width <= 3000){
        w = input$ref_inf_single_width
      }else{
        w = 3000
      }

      if(input$ref_inf_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_inf_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_inf_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      p = immune_all_datasets(datasets_rna_wes = datasets_rna_wes,dataset = ref_single_para()$datasource,gene = ref_single_para()$gene,immune_module = "immune_infiltration",dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
      print(p)
      dev.off()
      shinyjs::enable("ref_inf_single_down")
    }
  )

  output$immune_signature_datasets_all = renderPlotly({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    shinyjs::disable("sec")
    shinyjs::removeClass(id = "button_id",class = "fa fa-arrow-alt-circle-up")
    shinyjs::addClass(id = "button_id",class = "fa fa-spinner fa-spin")

    p = immune_all_datasets(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_pathway",dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)
    tmp = ggplotly(p) %>%
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0), height = 800) %>%
            rangeslider(start = 0,end = nrow(datasets_rna_wes[[ds]]$immune_pathway)+1) %>% plotly::config(displayModeBar = FALSE)

    shinyjs::enable("sec")
    shinyjs::removeClass(id = "button_id",class = "fa fa-spinner fa-spin")
    shinyjs::addClass(id = "button_id",class = "fa fa-arrow-alt-circle-up")

    return(tmp)
  }) %>% bindCache(ref_single_para()$gene,ref_single_para()$datasource,ref_cohort())

  output$ref_sig_single_down = downloadHandler(
    filename = function(){
      if(input$ref_sig_single_file == 'pdf'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","sig", '.',"pdf")
      }else if(input$ref_sig_single_file == 'png'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","sig", '.',"png")
      }else if(input$ref_sig_single_file == 'jpeg'){
        paste0(ref_single_para()$gene,"_",input$dataset,".","sig", '.',"jpeg")
      }else{
        paste0(ref_single_para()$gene,"_",input$dataset,".","sig", '.',"tiff")
      }

    },
    content = function(file){
      shinyjs::disable("ref_sig_single_down")
      if(input$ref_sig_single_res <= 300){
        r = input$ref_sig_single_res
      }else{
        r = 300
      }
      if(input$ref_sig_single_file == "pdf" & input$ref_sig_single_height <= 50){
        h = input$ref_sig_single_height
      }else if(input$ref_sig_single_file == "pdf" & input$ref_sig_single_height > 50){
        h = 50
      }else if(input$ref_sig_single_file != "pdf" & input$ref_sig_single_height <= 3000){
        h = input$ref_sig_single_height
      }else{
        h = 3000
      }
      if(input$ref_sig_single_file == "pdf" & input$ref_sig_single_width <= 50){
        w = input$ref_sig_single_width
      }else if(input$ref_sig_single_file == "pdf" & input$ref_sig_single_width > 50){
        w = 50
      }else if(input$ref_sig_single_file != "pdf" & input$ref_sig_single_width <= 3000){
        w = input$ref_sig_single_width
      }else{
        w = 3000
      }

      if(input$ref_sig_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$ref_sig_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$ref_sig_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      p = immune_all_datasets(datasets_rna_wes = datasets_rna_wes,dataset = ref_single_para()$datasource,gene = ref_single_para()$gene,immune_module = "immune_pathway",dataset_mu = datasets_mu[[ref_single_para()$datasource]],mut = ref_cohort()$mut,wt = ref_cohort()$wt)
      print(p)
      dev.off()
      shinyjs::enable("ref_sig_single_down")
    }
  )

  output$immune_infiltration_datasets_one = renderPlot({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),FALSE))

    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");

    tmp = immune_one_datasets(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_infiltration",selection = input$immune_infiltration_datasets_input,outlier = input$outlier_datasets,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)

    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");

    return(tmp)
    }) %>% bindCache(ref_single_para()$gene,ref_single_para()$datasource,ref_cohort(),input$immune_infiltration_datasets_input,input$outlier_datasets)

  output$immune_signature_datasets_one = renderPlot({
    gene = ref_single_para()$gene
    ds = ref_single_para()$datasource
    ref_cohort = ref_cohort()
    validate(need(length(ref_cohort$mut) >=3 & length(ref_cohort$wt) >= 3,message = FALSE))
    validate(need(ds %in% c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),message = FALSE))

    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");

    tmp = immune_one_datasets(datasets_rna_wes = datasets_rna_wes,dataset = ds,gene = gene,immune_module = "immune_pathway",selection = input$immune_signature_datasets_input,outlier = input$outlier_datasets,dataset_mu = datasets_mu[[ds]],mut = ref_cohort$mut,wt = ref_cohort$wt)

    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");

    return(tmp)
  }) %>% bindCache(ref_single_para()$gene,ref_single_para()$datasource,ref_cohort(),input$immune_signature_datasets_input,input$outlier_datasets)
}