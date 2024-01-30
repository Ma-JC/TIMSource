TCGA_subtype_server <- function(input,output,session,cancer_type_subtype,Mutation_subtype_tcga,gene_tcga_subtype,TCGA,pathway_database,Mut_type_tcga_subtype,Wild_type_tcga_subtype,MT_OR_WT_tcga,min.pct_subtype,FC_subtype,pvalue_subtype){
  
  
  
  TCGA_cohort_subtype_before = reactive({
    Mutation_subtype = Mutation_subtype_tcga()
    cancer_type_subtype = cancer_type_subtype()
    Mut_type_tcga_subtype = Mut_type_tcga_subtype()
    Wild_type_tcga_subtype = Wild_type_tcga_subtype()
    
    TCGA_cohort_cal_subtype(TCGA = TCGA,Mutation_subtype = Mutation_subtype,cancer_type = cancer_type_subtype,Mut_type = Mut_type_tcga_subtype,Wild_type = Wild_type_tcga_subtype)
    
    })%>% bindCache(Mutation_subtype_tcga(),cancer_type_subtype(),Mut_type_tcga_subtype(),Wild_type_tcga_subtype())
  
  TCGA_cohort_subtype = reactive({
    TCGA_cohort_subtype_before = TCGA_cohort_subtype_before()

      validate(need(length(TCGA_cohort_subtype_before$mut) >=6 & length(TCGA_cohort_subtype_before$wt) >= 6,message = FALSE))
    
    TCGA_cohort_subtype_before
    
  }) %>% bindCache(TCGA_cohort_subtype_before())
  
  #############immune_infiltration###############
  
  output$inf_tcga_subtype_uidown = renderUI({
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    tagList(
      div(selectInput(inputId = "TCGA_inf_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_inf_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_inf_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_inf_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_inf_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_inf_subtype_file,{
    if(input$TCGA_inf_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_inf_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_inf_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_inf_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_inf_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_inf_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_inf_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec2_subtype,{
    TCGA_cohort_subtype_before = TCGA_cohort_subtype_before()
    
    if(
      (length(TCGA_cohort_subtype_before$mut) < 6) | (length(TCGA_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning_tcga_subtype_id")
      createAlert(session, "warning_tcga_subtype", "warning_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning_tcga_subtype_id")
    }
  })
  
  
  output$immune_infiltration1_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_tcga()
    genes = gene_tcga_subtype()
    cancer_type_subtype = cancer_type_subtype()
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_infiltration_subtype_tmpdata_compare(TCGA = TCGA,cancer_type = cancer_type_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
    
    
    if(table(tmp_data$module,tmp_data$immune_infiltration)[2,1] <6 | table(tmp_data$module,tmp_data$immune_infiltration)[3,1] < 6){
      closeAlert(session,"warning_tcga_subtype_id")
      createAlert(session, "warning_tcga_subtype", "warning_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning_tcga_subtype_id")
    }
    validate(need(table(tmp_data$module,tmp_data$immune_infiltration)[2,1] >=6 | table(tmp_data$module,tmp_data$immune_infiltration)[3,1] >=6,message = FALSE))
    
    

    
    p = immune_infiltration_subtype_compare(tmp_data_compare = tmp_data,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(p)
    
    
  }) %>% bindCache(Mutation_subtype_tcga(),gene_tcga_subtype(),cancer_type_subtype(),TCGA_cohort_subtype())
  
  
  output$TCGA_inf_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_tcga_subtype()
      Mutation_subtype = Mutation_subtype_tcga()
      
      if(input$TCGA_inf_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","inf", '.',"pdf")
      }else if(input$TCGA_inf_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","inf", '.',"png")
      }else if(input$TCGA_inf_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","inf", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","inf", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_tcga()
      genes = gene_tcga_subtype()
      cancer_type_subtype = cancer_type_subtype()
      TCGA_cohort_subtype = TCGA_cohort_subtype()
      
      if(input$TCGA_inf_subtype_res <= 300){
        r = input$TCGA_inf_subtype_res
      }else{
        r = 300
      }
      if(input$TCGA_inf_subtype_file == "pdf" & input$TCGA_inf_subtype_height <= 50){
        h = input$TCGA_inf_subtype_height
      }else if(input$TCGA_inf_subtype_file == "pdf" & input$TCGA_inf_subtype_height > 50){
        h = 50
      }else if(input$TCGA_inf_subtype_file != "pdf" & input$TCGA_inf_subtype_height <= 3000){
        h = input$TCGA_inf_subtype_height
      }else{
        h = 3000
      }
      if(input$TCGA_inf_subtype_file == "pdf" & input$TCGA_inf_subtype_width <= 50){
        w = input$TCGA_inf_subtype_width
      }else if(input$TCGA_inf_subtype_file == "pdf" & input$TCGA_inf_subtype_width > 50){
        w = 50
      }else if(input$TCGA_inf_subtype_file != "pdf" & input$TCGA_inf_subtype_width <= 3000){
        w = input$TCGA_inf_subtype_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_inf_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_inf_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_inf_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tmp_data = immune_infiltration_subtype_tmpdata_compare(TCGA = TCGA,cancer_type = cancer_type_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
      print(immune_infiltration_subtype_compare(tmp_data_compare = tmp_data,genes = genes))
      
      dev.off()
    }
  )
  
  output$immune_infiltration2_subtype = renderPlotly({
    Mutation_subtype = Mutation_subtype_tcga()
    genes = gene_tcga_subtype()
    cancer_type_subtype = cancer_type_subtype()
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    
    tmp_data = immune_one_subtype_tmpdata_compare(TCGA = TCGA,Mutation_subtype = Mutation_subtype,genes = genes,cancer_type = cancer_type_subtype,immune_module = "immune_infiltration",selection = input$immune_cell_type_subtype,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
    
    
    validate(need(table(tmp_data$module)[2] >=6 | table(tmp_data$module)[3] >=6,message = FALSE))
    
    
    p = immune_one_subtype_compare(tmp_data_compare = tmp_data,selection = input$immune_cell_type_subtype,genes = genes)
    tmp = ggplotly(p) %>% 
      layout(boxmode = "group",legend = list(orientation = "v",x =1.05,y = 0.5,xanchor = "center"))

    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(Mutation_subtype_tcga(),gene_tcga_subtype(),cancer_type_subtype(),TCGA_cohort_subtype(),input$immune_cell_type_subtype)
  
  #############immune_signature##################
  
  output$sig_tcga_subtype_uidown = renderUI({
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    tagList(
      div(selectInput(inputId = "TCGA_sig_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_sig_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sig_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sig_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_sig_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_sig_subtype_file,{
    if(input$TCGA_sig_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_sig_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_sig_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_sig_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_sig_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sig_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sig_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec2_subtype,{
    TCGA_cohort_subtype_before = TCGA_cohort_subtype_before()

    if(
      (length(TCGA_cohort_subtype_before$mut) < 6) | (length(TCGA_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning2_tcga_subtype_id")
      createAlert(session, "warning2_tcga_subtype", "warning2_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning2_tcga_subtype_id")
    }
  })
  
  
  output$immune_signature1_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_tcga()
    genes = gene_tcga_subtype()
    cancer_type_subtype = cancer_type_subtype()
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_signature_subtype_tmpdata_compare(TCGA = TCGA,cancer_type = cancer_type_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
    
    
    if(table(tmp_data$module,tmp_data$immune_pathway)[2,1] <6 | table(tmp_data$module,tmp_data$immune_pathway)[3,1] < 6){
      closeAlert(session,"warning2_tcga_subtype_id")
      createAlert(session, "warning2_tcga_subtype", "warning2_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning2_tcga_subtype_id")
    }
    validate(need(table(tmp_data$module,tmp_data$immune_pathway)[2,1] >=6 & table(tmp_data$module,tmp_data$immune_pathway)[3,1] >=6 ,message = FALSE))
    
    
    p = immune_signature_subtype_compare(tmp_data_compare = tmp_data,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(p)
    
    
  }) %>% bindCache(Mutation_subtype_tcga(),gene_tcga_subtype(),cancer_type_subtype(),TCGA_cohort_subtype())
  
  output$TCGA_sig_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_tcga_subtype()
      Mutation_subtype = Mutation_subtype_tcga()
      
      if(input$TCGA_sig_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sig", '.',"pdf")
      }else if(input$TCGA_sig_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sig", '.',"png")
      }else if(input$TCGA_sig_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sig", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sig", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_tcga()
      genes = gene_tcga_subtype()
      cancer_type_subtype = cancer_type_subtype()
      TCGA_cohort_subtype = TCGA_cohort_subtype()
      
      if(input$TCGA_sig_subtype_res <= 300){
        r = input$TCGA_sig_subtype_res
      }else{
        r = 300
      }
      if(input$TCGA_sig_subtype_file == "pdf" & input$TCGA_sig_subtype_height <= 50){
        h = input$TCGA_sig_subtype_height
      }else if(input$TCGA_sig_subtype_file == "pdf" & input$TCGA_sig_subtype_height > 50){
        h = 50
      }else if(input$TCGA_sig_subtype_file != "pdf" & input$TCGA_sig_subtype_height <= 3000){
        h = input$TCGA_sig_subtype_height
      }else{
        h = 3000
      }
      if(input$TCGA_sig_subtype_file == "pdf" & input$TCGA_sig_subtype_width <= 50){
        w = input$TCGA_sig_subtype_width
      }else if(input$TCGA_sig_subtype_file == "pdf" & input$TCGA_sig_subtype_width > 50){
        w = 50
      }else if(input$TCGA_sig_subtype_file != "pdf" & input$TCGA_sig_subtype_width <= 3000){
        w = input$TCGA_sig_subtype_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_sig_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_sig_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_sig_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      tmp_data = immune_signature_subtype_tmpdata_compare(TCGA = TCGA,cancer_type = cancer_type_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
      print(immune_signature_subtype_compare(tmp_data_compare = tmp_data,genes = genes))      
      dev.off()
    }
  )
  
 
  output$immune_signature2_subtype = renderPlotly({
    Mutation_subtype = Mutation_subtype_tcga()
    genes = gene_tcga_subtype()
    cancer_type_subtype = cancer_type_subtype()
    TCGA_cohort_subtype = TCGA_cohort_subtype()

    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_one_subtype_tmpdata_compare(TCGA = TCGA,Mutation_subtype = Mutation_subtype,genes = genes,cancer_type = cancer_type_subtype,immune_module = "immune_pathway",selection = input$immune_signature_subtype,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
    validate(need(table(tmp_data$module)[2] >=6 | table(tmp_data$module)[3] >=6,message = FALSE))
    
    
    p = immune_one_subtype_compare(tmp_data_compare = tmp_data,selection = input$immune_signature_subtype,genes = genes)
    tmp = ggplotly(p) %>% 
      layout(boxmode = "group",legend = list(orientation = "v",x =1.05,y = 0.5,xanchor = "center"))
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  })  %>% bindCache(Mutation_subtype_tcga(),gene_tcga_subtype(),cancer_type_subtype(),TCGA_cohort_subtype(),input$immune_signature_subtype)
  
  ############## DEG ###################
  
  output$DEG_tcga_subtype_uitabledown = renderUI({
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    downloadButton('TCGA_diff_subtype_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec2_subtype,{
    TCGA_cohort_subtype_before = TCGA_cohort_subtype_before()
    
    if(
      (length(TCGA_cohort_subtype_before$mut) < 6) | (length(TCGA_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning3_tcga_subtype_id")
      createAlert(session, "warning3_tcga_subtype", "warning3_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning3_tcga_subtype_id")
    }
  })
  
  
  DEG_table_subtype = reactive({
    Mutation_subtype = Mutation_subtype_tcga()
    genes = gene_tcga_subtype()
    cancer_type_subtype = cancer_type_subtype()
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    
    MT_OR_WT = MT_OR_WT_tcga()
    min.pct_subtype = min.pct_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = DEG_subtype_tmpdata(TCGA = TCGA,MT_OR_WT = MT_OR_WT,genes = genes,cancer_type = cancer_type_subtype,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
    
    
    
    if(nrow(tmp_data$tmp_data) < 6){
      closeAlert(session,"warning3_tcga_subtype_id")
      createAlert(session, "warning3_tcga_subtype", "warning3_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)

      closeAlert(session,"warning4_tcga_subtype_id")
      createAlert(session, "warning4_tcga_subtype", "warning4_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)

    }else{
      closeAlert(session,"warning3_tcga_subtype_id")
      closeAlert(session,"warning4_tcga_subtype_id")
    }
    validate(need(nrow(tmp_data$tmp_data) >=6,message = FALSE))
    
    
    tmp = DEG_subtype(TCGA = TCGA,min.pct = min.pct_subtype,tmp_data = tmp_data)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(Mutation_subtype_tcga(),gene_tcga_subtype(),cancer_type_subtype(),TCGA_cohort_subtype(),MT_OR_WT_tcga(),min.pct_subtype())
  
  output$DEG_tab_subtype = renderReactable({
    DEG_table_subtype = DEG_table_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_subtype,3),
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
  }) %>% bindCache(DEG_table_subtype())
  
  
  output$TCGA_diff_subtype_tabdown = downloadHandler(
    
    filename = function(){
      gene = Mutation_subtype_tcga()
      paste0(gene,"_",input$Cancer_type_subtype,".","diff", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = DEG_table_subtype(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  output$volcano_subtype = renderPlotly({
    
    tab = DEG_table_subtype()
    FC_subtype = FC_subtype()
    pvalue_subtype = pvalue_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_subtype,pvalue = pvalue_subtype)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(DEG_table_subtype(),FC_subtype(),pvalue_subtype())
  
  #############TCGA GSEA##############
  
  output$GSEA_tcga_subtype_uitabledown = renderUI({
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    downloadButton('TCGA_gsea_subtype_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec2_subtype,{
    TCGA_cohort_subtype_before = TCGA_cohort_subtype_before()

    if(
      (length(TCGA_cohort_subtype_before$mut) < 6) | (length(TCGA_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning4_tcga_subtype_id")
      createAlert(session, "warning4_tcga_subtype", "warning4_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning4_tcga_subtype_id")
    }
  })
  
  status_file <- tempfile()
  write("", status_file)
  
  tmp_gsea_tcga_subtype = reactive({
    diff = DEG_table_subtype()
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
    
    installr::kill_pid(pid = scan(status_file))
    tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_tcga_subtype]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
    write(tmp$get_pid(), status_file)
    
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
  output$GSEA_tab_tcga_subtype = renderReactable({
    
    if(tmp_gsea_tcga_subtype()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader() == 0){
        
        waiter_show(id = "tcga_subtype_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader(one_loader()+1)
        return(NULL)
        
        
      }else if(one_loader() == 1){
        
        shinyjs::disable(id = "TCGA_gsea_subtype_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "tcga_subtype_GSEA_loader")
      
      tmp_gsea_tcga_subtype = tmp_gsea_tcga_subtype()$get_result()@result
      tmp_gsea_tcga_subtype$ID = as.factor(tmp_gsea_tcga_subtype$ID)
      tmp_gsea_tcga_subtype$Description = as.factor(tmp_gsea_tcga_subtype$Description)
      tmp_gsea_tcga_subtype$Plot = NA
      tmp_gsea_tcga_subtype = tmp_gsea_tcga_subtype[,c(ncol(tmp_gsea_tcga_subtype),3:(ncol(tmp_gsea_tcga_subtype)-2))]
      tmp = reactable( tmp_gsea_tcga_subtype,
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
                           Shiny.setInputValue('tcga_subtype_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "TCGA_gsea_subtype_tabdown")
      one_loader(0)
      write("", status_file)
      return(tmp)
    }

    
  })
  
  observeEvent(input$tcga_subtype_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "TCGA_gsea_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "TCGA_gsea_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "TCGA_gsea_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "TCGA_gsea_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "TCGA_gsea_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_tcga_subtype"),type = 1)
    ))
  })
  
  observeEvent(input$TCGA_gsea_subtype_file,{
    if(input$TCGA_gsea_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_gsea_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_gsea_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_gsea_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_gsea_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_gsea_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_gsea_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$TCGA_gsea_subtype_tabdown = downloadHandler(
    
    filename = function(){
      gene = Mutation_subtype_tcga()
      paste0(gene,"_",input$Cancer_type_subtype,".","gsea", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = tmp_gsea_tcga_subtype()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  row_num_tcga_subtype = eventReactive(input$tcga_subtype_GSEA_show,{input$tcga_subtype_GSEA_show$index})
  
  output$GSEA_plot_tcga_subtype = renderPlot({
    
    my_gseaplot2(tmp_gsea_tcga_subtype()$get_result(),geneSetID = tmp_gsea_tcga_subtype()$get_result()@result$ID[row_num_tcga_subtype()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_tcga_subtype()$get_result()@result$ID[row_num_tcga_subtype()])
  })
  
  output$TCGA_gsea_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_tcga_subtype()
      Mutation_subtype = Mutation_subtype_tcga()
      MT_OR_WT_tcga = MT_OR_WT_tcga()
      if(input$TCGA_gsea_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT_tcga,")",".","gsea", '.',"pdf")
      }else if(input$TCGA_gsea_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT_tcga,")",".","gsea", '.',"png")
      }else if(input$TCGA_gsea_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT_tcga,")",".","gsea", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT_tcga,")",".","gsea", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_tcga()
      genes = gene_tcga_subtype()
      cancer_type_subtype = cancer_type_subtype()
      TCGA_cohort_subtype = TCGA_cohort_subtype()
      
      if(input$TCGA_gsea_subtype_res <= 300){
        r = input$TCGA_gsea_subtype_res
      }else{
        r = 300
      }
      if(input$TCGA_gsea_subtype_file == "pdf" & input$TCGA_gsea_subtype_height <= 50){
        h = input$TCGA_gsea_subtype_height
      }else if(input$TCGA_gsea_subtype_file == "pdf" & input$TCGA_gsea_subtype_height > 50){
        h = 50
      }else if(input$TCGA_gsea_subtype_file != "pdf" & input$TCGA_gsea_subtype_height <= 3000){
        h = input$TCGA_gsea_subtype_height
      }else{
        h = 3000
      }
      if(input$TCGA_gsea_subtype_file == "pdf" & input$TCGA_gsea_subtype_width <= 50){
        w = input$TCGA_gsea_subtype_width
      }else if(input$TCGA_gsea_subtype_file == "pdf" & input$TCGA_gsea_subtype_width > 50){
        w = 50
      }else if(input$TCGA_gsea_subtype_file != "pdf" & input$TCGA_gsea_subtype_width <= 3000){
        w = input$TCGA_gsea_subtype_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_gsea_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_gsea_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_gsea_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(
        
        my_gseaplot2(tmp_gsea_tcga_subtype()$get_result(),geneSetID = tmp_gsea_tcga_subtype()$get_result()@result$ID[row_num_tcga_subtype()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_tcga_subtype()$get_result()@result$ID[row_num_tcga_subtype()])
        
      )
      
      dev.off()
    }
  )
  
  #############Survival##################
  
  output$sur_tcga_subtype_uidown = renderUI({
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    tagList(
      div(selectInput(inputId = "TCGA_sur_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_sur_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sur_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sur_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_sur_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_sur_subtype_file,{
    if(input$TCGA_sur_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_sur_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_sur_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_sur_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_sur_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sur_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sur_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec2_subtype,{
    TCGA_cohort_subtype_before = TCGA_cohort_subtype_before()
    
    if(
      (length(TCGA_cohort_subtype_before$mut) < 6) | (length(TCGA_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning5_tcga_subtype_id")
      createAlert(session, "warning5_tcga_subtype", "warning5_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning5_tcga_subtype_id")
    }
  })
  
  
  
  height = reactive({    if(cancer_type_subtype() == "LAML"){350}else{1100}  }) %>% bindCache(cancer_type_subtype())
  
  output$TCGA_survival_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_tcga()
    genes = gene_tcga_subtype()
    cancer_type_subtype = cancer_type_subtype()
    TCGA_cohort_subtype = TCGA_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    
    tmp_data = TCGA_survival_subtype_tmpdata_compare(TCGA = TCGA,cancer_type = cancer_type_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
    
    
    
    if(table(tmp_data$module)[2] < 6 | table(tmp_data$module)[3] < 6){
      closeAlert(session,"warning5_tcga_subtype_id")
      createAlert(session, "warning5_tcga_subtype", "warning5_tcga_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning5_tcga_subtype_id")
    }
    validate(need(table(tmp_data$module)[2] >= 6 | table(tmp_data$module)[3] >= 6,message = FALSE))
    
    tmp = TCGA_survival_subtype_compare(tmp_data_compare = tmp_data,cancer_type = cancer_type_subtype,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  },width = 1100,height = height) %>% bindCache(Mutation_subtype_tcga(),gene_tcga_subtype(),cancer_type_subtype(),TCGA_cohort_subtype())
  
  output$TCGA_sur_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_tcga_subtype()
      Mutation_subtype = Mutation_subtype_tcga()
      
      if(input$TCGA_sur_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"pdf")
      }else if(input$TCGA_sur_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"png")
      }else if(input$TCGA_sur_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_tcga()
      genes = gene_tcga_subtype()
      cancer_type_subtype = cancer_type_subtype()
      TCGA_cohort_subtype = TCGA_cohort_subtype()
      
      if(input$TCGA_sur_subtype_res <= 300){
        r = input$TCGA_sur_subtype_res
      }else{
        r = 300
      }
      if(input$TCGA_sur_subtype_file == "pdf" & input$TCGA_sur_subtype_height <= 50){
        h = input$TCGA_sur_subtype_height
      }else if(input$TCGA_sur_subtype_file == "pdf" & input$TCGA_sur_subtype_height > 50){
        h = 50
      }else if(input$TCGA_sur_subtype_file != "pdf" & input$TCGA_sur_subtype_height <= 3000){
        h = input$TCGA_sur_subtype_height
      }else{
        h = 3000
      }
      if(input$TCGA_sur_subtype_file == "pdf" & input$TCGA_sur_subtype_width <= 50){
        w = input$TCGA_sur_subtype_width
      }else if(input$TCGA_sur_subtype_file == "pdf" & input$TCGA_sur_subtype_width > 50){
        w = 50
      }else if(input$TCGA_sur_subtype_file != "pdf" & input$TCGA_sur_subtype_width <= 3000){
        w = input$TCGA_sur_subtype_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_sur_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_sur_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_sur_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tmp_data = TCGA_survival_subtype_tmpdata_compare(TCGA = TCGA,cancer_type = cancer_type_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = TCGA_cohort_subtype$mut,wt_patient_id = TCGA_cohort_subtype$wt)
      print(TCGA_survival_subtype_compare(tmp_data_compare = tmp_data,cancer_type = cancer_type_subtype,genes = genes))
      
      dev.off()
    }
  )
  
}