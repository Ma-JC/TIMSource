CPTAC_subtype_server <- function(input,output,session,cancer_type2_subtype,Mutation_subtype_cptac,gene_cptac_subtype,CPTAC,pathway_database,Mut_type_cptac_subtype,Wild_type_cptac_subtype,
                                 MT_OR_WT_cptac_rna,min.pct_CPTAC_rna_subtype,FC_CPTAC_rna_subtype,pvalue_CPTAC_rna_subtype,MT_OR_WT_cptac_protein,min.pct_CPTAC_protein_subtype,FC_CPTAC_protein_subtype,pvalue_CPTAC_protein_subtype){
  
  CPTAC_cohort_subtype_before = reactive({
    Mutation_subtype = Mutation_subtype_cptac()
    cancer_type2_subtype = cancer_type2_subtype()
    Mut_type_cptac_subtype = Mut_type_cptac_subtype()
    Wild_type_cptac_subtype = Wild_type_cptac_subtype()
    
    CPTAC_cohort_cal_subtype(TCGA = CPTAC,Mutation_subtype = Mutation_subtype,cancer_type = cancer_type2_subtype,Mut_type = Mut_type_cptac_subtype,Wild_type = Wild_type_cptac_subtype)
  }) %>% bindCache(Mutation_subtype_cptac(),cancer_type2_subtype(),Mut_type_cptac_subtype(),Wild_type_cptac_subtype())
  
  CPTAC_cohort_subtype = reactive({
    CPTAC_cohort_subtype_before = CPTAC_cohort_subtype_before()
    
    validate(need(length(CPTAC_cohort_subtype_before$mut) >=6 & length(CPTAC_cohort_subtype_before$wt) >= 6,message = FALSE))
    
    CPTAC_cohort_subtype_before
    
  }) %>% bindCache(CPTAC_cohort_subtype_before())
   
  ################################immune_infiltration#################################
  
  output$infr_cptac_subtype_uidown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    tagList(
      div(selectInput(inputId = "CPTAC_infr_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_infr_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infr_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infr_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_infr_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_infr_subtype_file,{
    if(input$CPTAC_infr_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_infr_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_infr_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_infr_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_infr_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infr_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infr_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$infp_cptac_subtype_uidown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    validate(need(length(CPTAC_cohort_subtype$mut) >=3 & length(CPTAC_cohort_subtype$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_infp_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_infp_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infp_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infp_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_infp_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_infp_subtype_file,{
    if(input$CPTAC_infp_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_infp_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_infp_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_infp_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_infp_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infp_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infp_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec3_subtype,{
    CPTAC_cohort_subtype_before = CPTAC_cohort_subtype_before()
    
    if(
      (length(CPTAC_cohort_subtype_before$mut) < 6) | (length(CPTAC_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning_cptac_subtype_id")
      createAlert(session, "warning_cptac_subtype", "warning_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning_cptac_subtype_id")
    }
  })
  
  
  output$immune_infiltration1_CPTAC_rna_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_infiltration_CPTAC_rna_subtype_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    if(table(tmp_data$module,tmp_data$immune_infiltration)[2,1] <6 | table(tmp_data$module,tmp_data$immune_infiltration)[3,1] < 6){
      closeAlert(session,"warning_cptac_subtype_id")
      createAlert(session, "warning_cptac_subtype", "warning_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning_cptac_subtype_id")
    }
    validate(need(table(tmp_data$module,tmp_data$immune_infiltration)[2,1] >=6 | table(tmp_data$module,tmp_data$immune_infiltration)[3,1] >=6,message = FALSE))
    
    p = immune_infiltration_CPTAC_rna_subtype_compare(tmp_data_compare = tmp_data,genes = genes)

    
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(p)
    
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype())
  
  output$CPTAC_infr_subtype_down = downloadHandler(
    filename = function(){
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      
      if(input$CPTAC_infr_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infr", '.',"pdf")
      }else if(input$CPTAC_infr_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infr", '.',"png")
      }else if(input$CPTAC_infr_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infr", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infr", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_cptac()
      genes = gene_cptac_subtype()
      cancer_type2_subtype = cancer_type2_subtype()
      CPTAC_cohort_subtype = CPTAC_cohort_subtype()
      
      if(input$CPTAC_infr_subtype_res <= 300){
        r = input$CPTAC_infr_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_infr_subtype_file == "pdf" & input$CPTAC_infr_subtype_height <= 50){
        h = input$CPTAC_infr_subtype_height
      }else if(input$CPTAC_infr_subtype_file == "pdf" & input$CPTAC_infr_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_infr_subtype_file != "pdf" & input$CPTAC_infr_subtype_height <= 3000){
        h = input$CPTAC_infr_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_infr_subtype_file == "pdf" & input$CPTAC_infr_subtype_width <= 50){
        w = input$CPTAC_infr_subtype_width
      }else if(input$CPTAC_infr_subtype_file == "pdf" & input$CPTAC_infr_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_infr_subtype_file != "pdf" & input$CPTAC_infr_subtype_width <= 3000){
        w = input$CPTAC_infr_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_infr_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_infr_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_infr_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tmp_data = immune_infiltration_CPTAC_rna_subtype_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
      print(immune_infiltration_CPTAC_rna_subtype_compare(tmp_data_compare = tmp_data,genes = genes))
      
      dev.off()
    }
  )
  
  
  output$immune_infiltration1_CPTAC_protein_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_infiltration_CPTAC_protein_subtype_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    if(table(tmp_data$module,tmp_data$immune_infiltration_protein)[2,1] <6 | table(tmp_data$module,tmp_data$immune_infiltration_protein)[3,1] < 6){
      closeAlert(session,"warning_cptac_subtype_id")
      createAlert(session, "warning_cptac_subtype", "warning_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning_cptac_subtype_id")
    }
    validate(need(table(tmp_data$module,tmp_data$immune_infiltration_protein)[2,1] >=6 | table(tmp_data$module,tmp_data$immune_infiltration_protein)[3,1] >=6,message = FALSE))

    
    
    p = immune_infiltration_CPTAC_protein_subtype_compare(tmp_data_compare = tmp_data,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(p)
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype())
  
  
  output$CPTAC_infp_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      
      if(input$CPTAC_infp_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infp", '.',"pdf")
      }else if(input$CPTAC_infp_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infp", '.',"png")
      }else if(input$CPTAC_infp_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infp", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","infp", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_cptac()
      genes = gene_cptac_subtype()
      cancer_type2_subtype = cancer_type2_subtype()
      CPTAC_cohort_subtype = CPTAC_cohort_subtype()
      
      if(input$CPTAC_infp_subtype_res <= 300){
        r = input$CPTAC_infp_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_infp_subtype_file == "pdf" & input$CPTAC_infp_subtype_height <= 50){
        h = input$CPTAC_infp_subtype_height
      }else if(input$CPTAC_infp_subtype_file == "pdf" & input$CPTAC_infp_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_infp_subtype_file != "pdf" & input$CPTAC_infp_subtype_height <= 3000){
        h = input$CPTAC_infp_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_infp_subtype_file == "pdf" & input$CPTAC_infp_subtype_width <= 50){
        w = input$CPTAC_infp_subtype_width
      }else if(input$CPTAC_infp_subtype_file == "pdf" & input$CPTAC_infp_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_infp_subtype_file != "pdf" & input$CPTAC_infp_subtype_width <= 3000){
        w = input$CPTAC_infp_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_infp_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_infp_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_infp_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tmp_data = immune_infiltration_CPTAC_protein_subtype_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,Mutation_subtype = Mutation_subtype,genes = genes,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
      print(immune_infiltration_CPTAC_protein_subtype_compare(tmp_data_compare = tmp_data,genes = genes))
      
      dev.off()
    }
  )
  
  output$immune_infiltration2_CPTAC_rna_subtype = renderPlotly({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_one_CPTAC_subtype_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,Mutation_subtype = Mutation_subtype,genes = genes,immune_module = "immune_infiltration",selection = input$immune_cell_type_CPTAC_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    validate(need(table(tmp_data$module)[2] >=6 | table(tmp_data$module)[3] >=6,message = FALSE))
    
    p = immune_one_CPTAC_subtype_compare(tmp_data_compare = tmp_data,immune_module = "immune_infiltration",selection = input$immune_cell_type_CPTAC_subtype,genes = genes)
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
    
    }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype(),input$immune_cell_type_CPTAC_subtype)
  
  output$immune_infiltration2_CPTAC_protein_subtype = renderPlotly({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_one_CPTAC_subtype_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,Mutation_subtype = Mutation_subtype,genes = genes,immune_module = "immune_infiltration_protein",selection = input$immune_cell_type_CPTAC_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    validate(need(table(tmp_data$module)[2] >=6 | table(tmp_data$module)[3] >=6,message = FALSE))
    
    p = immune_one_CPTAC_subtype_compare(tmp_data_compare = tmp_data,immune_module = "immune_infiltration_protein",selection = input$immune_cell_type_CPTAC_subtype,genes = genes)
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
    
    }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype(),input$immune_cell_type_CPTAC_subtype)
  
  ##########################immune_signature###########################
  
  output$sigr_cptac_subtype_uidown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    tagList(
      div(selectInput(inputId = "CPTAC_sigr_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sigr_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigr_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigr_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sigr_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sigr_subtype_file,{
    if(input$CPTAC_sigr_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sigr_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sigr_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$sigp_cptac_subtype_uidown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    tagList(
      div(selectInput(inputId = "CPTAC_sigp_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sigp_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigp_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigp_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sigp_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sigp_subtype_file,{
    if(input$CPTAC_sigp_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sigp_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sigp_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec3_subtype,{
    CPTAC_cohort_subtype_before = CPTAC_cohort_subtype_before()

    if(
      (length(CPTAC_cohort_subtype_before$mut) < 6) | (length(CPTAC_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning2_cptac_subtype_id")
      createAlert(session, "warning2_cptac_subtype", "warning2_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning2_cptac_subtype_id")
    }
  })
  
  output$immune_signature1_CPTAC_rna_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_signature_CPTAC_rna_subtype_tmpdata_compare(TCGA = CPTAC,Mutation_subtype = Mutation_subtype,genes = genes,cancer_type = cancer_type2_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    if(table(tmp_data$module,tmp_data$immune_pathway)[2,1] <6 | table(tmp_data$module,tmp_data$immune_pathway)[3,1] < 6){
      closeAlert(session,"warning2_cptac_subtype_id")
      createAlert(session, "warning2_cptac_subtype", "warning2_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning2_cptac_subtype_id")
    }
    validate(need(table(tmp_data$module,tmp_data$immune_pathway)[2,1] >=6 | table(tmp_data$module,tmp_data$immune_pathway)[3,1] >=6,message = FALSE))
    
    p = immune_signature_CPTAC_rna_subtype_compare(tmp_data_compare = tmp_data,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(p)
    
    
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype())
  
  
  output$CPTAC_sigr_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      
      if(input$CPTAC_sigr_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigr", '.',"pdf")
      }else if(input$CPTAC_sigr_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigr", '.',"png")
      }else if(input$CPTAC_sigr_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigr", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigr", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_cptac()
      genes = gene_cptac_subtype()
      cancer_type2_subtype = cancer_type2_subtype()
      CPTAC_cohort_subtype = CPTAC_cohort_subtype()
      
      if(input$CPTAC_sigr_subtype_res <= 300){
        r = input$CPTAC_sigr_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_sigr_subtype_file == "pdf" & input$CPTAC_sigr_subtype_height <= 50){
        h = input$CPTAC_sigr_subtype_height
      }else if(input$CPTAC_sigr_subtype_file == "pdf" & input$CPTAC_sigr_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_sigr_subtype_file != "pdf" & input$CPTAC_sigr_subtype_height <= 3000){
        h = input$CPTAC_sigr_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sigr_subtype_file == "pdf" & input$CPTAC_sigr_subtype_width <= 50){
        w = input$CPTAC_sigr_subtype_width
      }else if(input$CPTAC_sigr_subtype_file == "pdf" & input$CPTAC_sigr_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_sigr_subtype_file != "pdf" & input$CPTAC_sigr_subtype_width <= 3000){
        w = input$CPTAC_sigr_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sigr_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sigr_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sigr_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tmp_data = immune_signature_CPTAC_rna_subtype_tmpdata_compare(TCGA = CPTAC,Mutation_subtype = Mutation_subtype,genes = genes,cancer_type = cancer_type2_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
      print(immune_signature_CPTAC_rna_subtype_compare(tmp_data_compare = tmp_data,genes = genes))
      
      dev.off()
    }
  )
  
  
  output$immune_signature1_CPTAC_protein_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_signature_CPTAC_protein_subtype_tmpdata_compare(TCGA = CPTAC,genes = genes,Mutation_subtype = Mutation_subtype,cancer_type = cancer_type2_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    if(table(tmp_data$module,tmp_data$immune_pathway_protein)[2,1] <6 | table(tmp_data$module,tmp_data$immune_pathway_protein)[3,1] < 6){
      closeAlert(session,"warning2_cptac_subtype_id")
      createAlert(session, "warning2_cptac_subtype", "warning2_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning2_cptac_subtype_id")
    }
    validate(need(table(tmp_data$module,tmp_data$immune_pathway_protein)[2,1] >=6 & table(tmp_data$module,tmp_data$immune_pathway_protein)[3,1] >= 6,message = FALSE))

    p = immune_signature_CPTAC_protein_subtype_compare(tmp_data_compare = tmp_data,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(p)
    
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype())
  
  
  output$CPTAC_sigp_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      
      if(input$CPTAC_sigp_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigp", '.',"pdf")
      }else if(input$CPTAC_sigp_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigp", '.',"png")
      }else if(input$CPTAC_sigp_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigp", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sigp", '.',"tiff")
      }
      
    },
    content = function(file){
      Mutation_subtype = Mutation_subtype_cptac()
      genes = gene_cptac_subtype()
      cancer_type2_subtype = cancer_type2_subtype()
      CPTAC_cohort_subtype = CPTAC_cohort_subtype()
      
      if(input$CPTAC_sigp_subtype_res <= 300){
        r = input$CPTAC_sigp_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_sigp_subtype_file == "pdf" & input$CPTAC_sigp_subtype_height <= 50){
        h = input$CPTAC_sigp_subtype_height
      }else if(input$CPTAC_sigp_subtype_file == "pdf" & input$CPTAC_sigp_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_sigp_subtype_file != "pdf" & input$CPTAC_sigp_subtype_height <= 3000){
        h = input$CPTAC_sigp_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sigp_subtype_file == "pdf" & input$CPTAC_sigp_subtype_width <= 50){
        w = input$CPTAC_sigp_subtype_width
      }else if(input$CPTAC_sigp_subtype_file == "pdf" & input$CPTAC_sigp_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_sigp_subtype_file != "pdf" & input$CPTAC_sigp_subtype_width <= 3000){
        w = input$CPTAC_sigp_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sigp_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sigp_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sigp_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tmp_data = immune_signature_CPTAC_protein_subtype_tmpdata_compare(TCGA = CPTAC,genes = genes,Mutation_subtype = Mutation_subtype,cancer_type = cancer_type2_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
      print(immune_signature_CPTAC_protein_subtype_compare(tmp_data_compare = tmp_data,genes = genes))
      
      dev.off()
    }
  )
  
  
  output$immune_signature2_CPTAC_rna_subtype = renderPlotly({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_one_CPTAC_subtype_tmpdata_compare(TCGA = CPTAC,Mutation_subtype = Mutation_subtype,genes = genes,cancer_type = cancer_type2_subtype,immune_module = "immune_pathway",selection = input$immune_signature_CPTAC_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    validate(need(table(tmp_data$module)[2] >=6 | table(tmp_data$module)[3] >=6,message = FALSE))
    
    p = immune_one_CPTAC_subtype_compare(tmp_data_compare = tmp_data,immune_module = "immune_pathway",selection = input$immune_signature_CPTAC_subtype,genes = genes)
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
    
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype(),input$immune_signature_CPTAC_subtype)
  
  output$immune_signature2_CPTAC_protein_subtype = renderPlotly({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    validate(need(input$immune_signature_CPTAC_subtype %in% rownames(CPTAC[[cancer_type2_subtype]][["immune_pathway_protein"]]),paste("The pathway",input$immune_signature_CPTAC_subtype,"does not exist in",cancer_type2_subtype,sep = " ")))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = immune_one_CPTAC_subtype_tmpdata_compare(TCGA = CPTAC,Mutation_subtype = Mutation_subtype,genes = genes,cancer_type = cancer_type2_subtype,immune_module = "immune_pathway_protein",selection = input$immune_signature_CPTAC_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    validate(need(table(tmp_data$module)[2] >=6 | table(tmp_data$module)[3] >=6,message = FALSE))
    
    
    p = immune_one_CPTAC_subtype_compare(tmp_data_compare = tmp_data,immune_module = "immune_pathway_protein",selection = input$immune_signature_CPTAC_subtype,genes = genes)
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
    
    })  %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype(),input$immune_signature_CPTAC_subtype)
  
  #####################################DEG######################################
  
  output$DEGr_cptac_subtype_uitabledown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    validate(need(length(CPTAC_cohort_subtype$mut) >=3 & length(CPTAC_cohort_subtype$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_diffr_subtype_tabdown',label = 'Download Table')
  })
  
  output$DEGp_cptac_subtype_uitabledown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    validate(need(length(CPTAC_cohort_subtype$mut) >=3 & length(CPTAC_cohort_subtype$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_diffp_subtype_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec3_subtype,{
    CPTAC_cohort_subtype_before = CPTAC_cohort_subtype_before()

    if(
      (length(CPTAC_cohort_subtype_before$mut) < 6) | (length(CPTAC_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning3_cptac_subtype_id")
      createAlert(session, "warning3_cptac_subtype", "warning3_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning3_cptac_subtype_id")
    }
  })
  
  
  DEG_table_CPTAC_rna_subtype = reactive({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    MT_OR_WT = MT_OR_WT_cptac_rna()
    min.pct_CPTAC_rna_subtype = min.pct_CPTAC_rna_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = DEG_CPTAC_rna_subtype_tmpdata(TCGA = CPTAC,MT_OR_WT = MT_OR_WT,genes = genes,cancer_type = cancer_type2_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    if(nrow(tmp_data$tmp_data) < 6 ){
      closeAlert(session,"warning3_cptac_subtype_id")
      createAlert(session, "warning3_cptac_subtype", "warning3_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
      
      closeAlert(session,"warning4_cptac_subtype_id")
      createAlert(session, "warning4_cptac_subtype", "warning4_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning3_cptac_subtype_id")
      closeAlert(session,"warning4_cptac_subtype_id")
    }
    validate(need(nrow(tmp_data$tmp_data) >=6,message = FALSE))
    
    tmp = DEG_CPTAC_rna_subtype(TCGA = CPTAC,min.pct = min.pct_CPTAC_rna_subtype,tmp_data = tmp_data)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype(),MT_OR_WT_cptac_rna(),min.pct_CPTAC_rna_subtype())
  
  output$DEG_tab_CPTAC_rna_subtype = renderReactable({
    DEG_table_CPTAC_rna_subtype = DEG_table_CPTAC_rna_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    
    tmp = reactable( round(DEG_table_CPTAC_rna_subtype,3),
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
  }) %>% bindCache(DEG_table_CPTAC_rna_subtype())
  
  output$CPTAC_diffr_subtype_tabdown = downloadHandler(
    
    filename = function(){
      gene = Mutation_subtype_cptac()
      paste0(gene,"_",input$Cancer_type2_subtype,".","diffr", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = DEG_table_CPTAC_rna_subtype(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  output$volcano_CPTAC_rna_subtype = renderPlotly({
    
    tab = DEG_table_CPTAC_rna_subtype()
    FC_CPTAC_rna_subtype = FC_CPTAC_rna_subtype()
    pvalue_CPTAC_rna_subtype = pvalue_CPTAC_rna_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_CPTAC_rna_subtype,pvalue = pvalue_CPTAC_rna_subtype)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(DEG_table_CPTAC_rna_subtype(),FC_CPTAC_rna_subtype(),pvalue_CPTAC_rna_subtype())
  
  DEG_table_CPTAC_protein_subtype = reactive({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    
    MT_OR_WT = MT_OR_WT_cptac_protein()
    min.pct_CPTAC_protein_subtype = min.pct_CPTAC_protein_subtype()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = DEG_CPTAC_protein_subtype_tmpdata(TCGA = CPTAC,MT_OR_WT = MT_OR_WT,genes = genes,cancer_type = cancer_type2_subtype,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    
    if(nrow(tmp_data$tmp_data) < 6 ){
      closeAlert(session,"warning3_cptac_subtype_id")
      createAlert(session, "warning3_cptac_subtype", "warning3_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
      
      closeAlert(session,"warning4_cptac_subtype_id")
      createAlert(session, "warning4_cptac_subtype", "warning4_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
      
    }else{
      closeAlert(session,"warning3_cptac_subtype_id")
      closeAlert(session,"warning4_cptac_subtype_id")
    }
    validate(need(nrow(tmp_data$tmp_data) >=6,message = FALSE))
    
    tmp = DEG_CPTAC_protein_subtype(TCGA = CPTAC,min.pct = min.pct_CPTAC_protein_subtype,tmp_data = tmp_data)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    
  }) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype(),MT_OR_WT_cptac_protein(),min.pct_CPTAC_protein_subtype())
  
  output$DEG_tab_CPTAC_protein_subtype = renderReactable({
    DEG_table_CPTAC_protein_subtype = DEG_table_CPTAC_protein_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_CPTAC_protein_subtype,3),
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
  }) %>% bindCache(DEG_table_CPTAC_protein_subtype())
  
  
  output$CPTAC_diffp_subtype_tabdown = downloadHandler(
    
    filename = function(){
      gene = Mutation_subtype_cptac()
      paste0(gene,"_",input$Cancer_type2_subtype,".","diffp", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = DEG_table_CPTAC_protein_subtype(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  output$volcano_CPTAC_protein_subtype = renderPlotly({
    
    tab = DEG_table_CPTAC_protein_subtype()
    FC_CPTAC_protein_subtype = FC_CPTAC_protein_subtype()
    pvalue_CPTAC_protein_subtype = pvalue_CPTAC_protein_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_CPTAC_protein_subtype,pvalue = pvalue_CPTAC_protein_subtype)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    
  }) %>% bindCache(DEG_table_CPTAC_protein_subtype(),FC_CPTAC_protein_subtype(),pvalue_CPTAC_protein_subtype())
  
  #################################################CPTAC GSEA######################################################
  
  output$GSEAr_cptac_subtype_uitabledown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    validate(need(length(CPTAC_cohort_subtype$mut) >=3 & length(CPTAC_cohort_subtype$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_gsear_subtype_tabdown',label = 'Download Table')
  })
  
  output$GSEAp_cptac_subtype_uitabledown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    validate(need(length(CPTAC_cohort_subtype$mut) >=3 & length(CPTAC_cohort_subtype$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_gseap_subtype_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec3_subtype,{
    CPTAC_cohort_subtype_before = CPTAC_cohort_subtype_before()

    if(
      (length(CPTAC_cohort_subtype_before$mut) < 6) | (length(CPTAC_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning4_cptac_subtype_id")
      createAlert(session, "warning4_cptac_subtype", "warning4_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning4_cptac_subtype_id")
    }
  })
  
  
  status_file_rna <- tempfile()
  write("", status_file_rna)
  status_file_protein <- tempfile()
  write("", status_file_protein)
  
  tmp_gsea_cptac_rna_subtype = reactive({
    diff = DEG_table_CPTAC_rna_subtype()
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
    
    installr::kill_pid(pid = scan(status_file_rna))
    tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_cptac_rna_subtype]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
    write(tmp$get_pid(), status_file_rna)
    
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
  
  one_loader_rna = reactiveVal(0)
  output$GSEA_tab_cptac_rna_subtype = renderReactable({
    
    if(tmp_gsea_cptac_rna_subtype()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader_rna() == 0){
        
        waiter_show(id = "cptac_subtype_rna_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader_rna(one_loader_rna()+1)
        return(NULL)
        
        
      }else if(one_loader_rna() == 1){
        
        shinyjs::disable(id = "CPTAC_gsear_subtype_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "cptac_subtype_rna_GSEA_loader")
      
      tmp_gsea_cptac_rna_subtype = tmp_gsea_cptac_rna_subtype()$get_result()@result
      tmp_gsea_cptac_rna_subtype$ID = as.factor(tmp_gsea_cptac_rna_subtype$ID)
      tmp_gsea_cptac_rna_subtype$Description = as.factor(tmp_gsea_cptac_rna_subtype$Description)
      tmp_gsea_cptac_rna_subtype$Plot = NA
      tmp_gsea_cptac_rna_subtype = tmp_gsea_cptac_rna_subtype[,c(ncol(tmp_gsea_cptac_rna_subtype),3:(ncol(tmp_gsea_cptac_rna_subtype)-2))]
      tmp = reactable( tmp_gsea_cptac_rna_subtype,
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
                           Shiny.setInputValue('cptac_subtype_rna_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "CPTAC_gsear_subtype_tabdown")
      one_loader_rna(0)
      write("", status_file_rna)
      return(tmp)
    }
    
    
    

  })
  
  observeEvent(input$cptac_subtype_rna_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "CPTAC_gsear_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "CPTAC_gsear_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gsear_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gsear_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "CPTAC_gsear_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_cptac_rna_subtype"),type = 1)
    ))
  })
  
  observeEvent(input$CPTAC_gsear_subtype_file,{
    if(input$CPTAC_gsear_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_gsear_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_gsear_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$CPTAC_gsear_subtype_tabdown = downloadHandler(
    
    filename = function(){
      gene = Mutation_subtype_cptac()
      paste0(gene,"_",input$Cancer_type2_subtype,".","gsear", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = tmp_gsea_cptac_rna_subtype()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  row_num_cptac_rna_subtype = eventReactive(input$cptac_subtype_rna_GSEA_show,{input$cptac_subtype_rna_GSEA_show$index})
  
  output$GSEA_plot_cptac_rna_subtype = renderPlot({
    
    
    my_gseaplot2(tmp_gsea_cptac_rna_subtype()$get_result(),geneSetID = tmp_gsea_cptac_rna_subtype()$get_result()@result$ID[row_num_cptac_rna_subtype()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_cptac_rna_subtype()$get_result()@result$ID[row_num_cptac_rna_subtype()])
    
    
  })
  
  
  output$CPTAC_gsear_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      MT_OR_WT = MT_OR_WT_cptac_rna()
      if(input$CPTAC_gsear_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gsear", '.',"pdf")
      }else if(input$CPTAC_gsear_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gsear", '.',"png")
      }else if(input$CPTAC_gsear_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gsear", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gsear", '.',"tiff")
      }
      
    },
    content = function(file){
      
      if(input$CPTAC_gsear_subtype_res <= 300){
        r = input$CPTAC_gsear_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_gsear_subtype_file == "pdf" & input$CPTAC_gsear_subtype_height <= 50){
        h = input$CPTAC_gsear_subtype_height
      }else if(input$CPTAC_gsear_subtype_file == "pdf" & input$CPTAC_gsear_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_gsear_subtype_file != "pdf" & input$CPTAC_gsear_subtype_height <= 3000){
        h = input$CPTAC_gsear_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_gsear_subtype_file == "pdf" & input$CPTAC_gsear_subtype_width <= 50){
        w = input$CPTAC_gsear_subtype_width
      }else if(input$CPTAC_gsear_subtype_file == "pdf" & input$CPTAC_gsear_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_gsear_subtype_file != "pdf" & input$CPTAC_gsear_subtype_width <= 3000){
        w = input$CPTAC_gsear_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_gsear_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_gsear_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_gsear_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      print(
        my_gseaplot2(tmp_gsea_cptac_rna_subtype()$get_result(),geneSetID = tmp_gsea_cptac_rna_subtype()$get_result()@result$ID[row_num_cptac_rna_subtype()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_cptac_rna_subtype()$get_result()@result$ID[row_num_cptac_rna_subtype()])
            )
      
      dev.off()
    }
  )
  
  
  tmp_gsea_cptac_protein_subtype = reactive({
    diff = DEG_table_CPTAC_protein_subtype()
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
    
    installr::kill_pid(pid = scan(status_file_protein))
    tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_cptac_protein_subtype]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
    write(tmp$get_pid(), status_file_protein)
    
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
  
  one_loader_protein = reactiveVal(0)
  output$GSEA_tab_cptac_protein_subtype = renderReactable({
    
    if(tmp_gsea_cptac_protein_subtype()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader_protein() == 0){
        
        waiter_show(id = "cptac_subtype_protein_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader_protein(one_loader_protein()+1)
        return(NULL)
        
        
      }else if(one_loader_protein() == 1){
        
        shinyjs::disable(id = "CPTAC_gseap_subtype_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "cptac_subtype_protein_GSEA_loader")
      
      tmp_gsea_cptac_protein_subtype = tmp_gsea_cptac_protein_subtype()$get_result()@result
      tmp_gsea_cptac_protein_subtype$ID = as.factor(tmp_gsea_cptac_protein_subtype$ID)
      tmp_gsea_cptac_protein_subtype$Description = as.factor(tmp_gsea_cptac_protein_subtype$Description)
      tmp_gsea_cptac_protein_subtype$Plot = NA
      tmp_gsea_cptac_protein_subtype = tmp_gsea_cptac_protein_subtype[,c(ncol(tmp_gsea_cptac_protein_subtype),3:(ncol(tmp_gsea_cptac_protein_subtype)-2))]
      tmp = reactable( tmp_gsea_cptac_protein_subtype,
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
                           Shiny.setInputValue('cptac_subtype_protein_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "CPTAC_gseap_subtype_tabdown")
      one_loader_protein(0)
      write("", status_file_protein)
      return(tmp)
    }

  })
  
  
  observeEvent(input$cptac_subtype_protein_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "CPTAC_gseap_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "CPTAC_gseap_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gseap_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gseap_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "CPTAC_gseap_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_cptac_protein_subtype"),type = 1)
    ))
  })
  
  observeEvent(input$CPTAC_gseap_subtype_file,{
    if(input$CPTAC_gseap_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_gseap_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_gseap_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$CPTAC_gseap_subtype_tabdown = downloadHandler(
    
    filename = function(){
      gene = Mutation_subtype_cptac()
      paste0(gene,"_",input$Cancer_type2_subtype,".","gseap", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = tmp_gsea_cptac_protein_subtype()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  row_num_cptac_protein_subtype = eventReactive(input$cptac_subtype_protein_GSEA_show,{input$cptac_subtype_protein_GSEA_show$index})
  
  output$GSEA_plot_cptac_protein_subtype = renderPlot({
    
    
    my_gseaplot2(tmp_gsea_cptac_protein_subtype()$get_result(),geneSetID = tmp_gsea_cptac_protein_subtype()$get_result()@result$ID[row_num_cptac_protein_subtype()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_cptac_protein_subtype()$get_result()@result$ID[row_num_cptac_protein_subtype()])
    
    
  })
  
  
  output$CPTAC_gseap_subtype_down = downloadHandler(
    
    filename = function(){
      
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      MT_OR_WT = MT_OR_WT_cptac_protein()
      if(input$CPTAC_gseap_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gseap", '.',"pdf")
      }else if(input$CPTAC_gseap_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gseap", '.',"png")
      }else if(input$CPTAC_gseap_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gseap", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype(",MT_OR_WT,")",".","gseap", '.',"tiff")
      }
      
    },
    content = function(file){
      
      if(input$CPTAC_gseap_subtype_res <= 300){
        r = input$CPTAC_gseap_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_gseap_subtype_file == "pdf" & input$CPTAC_gseap_subtype_height <= 50){
        h = input$CPTAC_gseap_subtype_height
      }else if(input$CPTAC_gseap_subtype_file == "pdf" & input$CPTAC_gseap_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_gseap_subtype_file != "pdf" & input$CPTAC_gseap_subtype_height <= 3000){
        h = input$CPTAC_gseap_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_gseap_subtype_file == "pdf" & input$CPTAC_gseap_subtype_width <= 50){
        w = input$CPTAC_gseap_subtype_width
      }else if(input$CPTAC_gseap_subtype_file == "pdf" & input$CPTAC_gseap_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_gseap_subtype_file != "pdf" & input$CPTAC_gseap_subtype_width <= 3000){
        w = input$CPTAC_gseap_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_gseap_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_gseap_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_gseap_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      print(
        
        my_gseaplot2(tmp_gsea_cptac_rna_subtype()$get_result(),geneSetID = tmp_gsea_cptac_rna_subtype()$get_result()@result$ID[row_num_cptac_rna_subtype()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_cptac_rna_subtype()$get_result()@result$ID[row_num_cptac_rna_subtype()])
        
      )
      
      dev.off()
    }
  )
  
  #############Survival##################
  
  output$sur_cptac_subtype_uidown = renderUI({
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    validate(need(length(CPTAC_cohort_subtype$mut) >=3 & length(CPTAC_cohort_subtype$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_sur_subtype_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sur_subtype_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sur_subtype_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sur_subtype_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sur_subtype_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sur_subtype_file,{
    if(input$CPTAC_sur_subtype_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sur_subtype_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sur_subtype_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sur_subtype_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sur_subtype_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sur_subtype_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sur_subtype_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec3_subtype,{
    CPTAC_cohort_subtype_before = CPTAC_cohort_subtype_before()

    if(
      (length(CPTAC_cohort_subtype_before$mut) < 6) | (length(CPTAC_cohort_subtype_before$wt) < 6)
    ){
      closeAlert(session,"warning5_cptac_subtype_id")
      createAlert(session, "warning5_cptac_subtype", "warning5_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning5_cptac_subtype_id")
    }
  })
  
  
  output$CPTAC_survival_subtype = renderPlot({
    Mutation_subtype = Mutation_subtype_cptac()
    genes = gene_cptac_subtype()
    cancer_type2_subtype = cancer_type2_subtype()
    CPTAC_cohort_subtype = CPTAC_cohort_subtype()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp_data = CPTAC_survival_subtype_fun_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,genes = genes,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
    if(nrow(tmp_data) < 6 ){
      closeAlert(session,"warning5_cptac_subtype_id")
      createAlert(session, "warning5_cptac_subtype", "warning5_cptac_subtype_id", title = "Warning",style = "danger",
                  content = "The number of patients with high or low expression of selected genes < 3", append = FALSE)
    }else{
      closeAlert(session,"warning5_cptac_subtype_id")
    }
    validate(need(nrow(tmp_data) >=6,message = FALSE))
    
    tmp = CPTAC_survival_subtype_fun_compare(tmp_data_compare = tmp_data,genes = genes)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  },width = 1100,height = 800) %>% bindCache(Mutation_subtype_cptac(),gene_cptac_subtype(),cancer_type2_subtype(),CPTAC_cohort_subtype())
  
  
  output$CPTAC_sur_subtype_down = downloadHandler(
    
    filename = function(){
      genes = gene_cptac_subtype()
      Mutation_subtype = Mutation_subtype_cptac()
      
      if(input$CPTAC_sur_subtype_file == 'pdf'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"pdf")
      }else if(input$CPTAC_sur_subtype_file == 'png'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"png")
      }else if(input$CPTAC_sur_subtype_file == 'jpeg'){
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"jpeg")
      }else{
        paste0(genes,"_",input$Cancer_type2_subtype,"_",Mutation_subtype,"_subtype",".","sur", '.',"tiff")
      }
      
    },
    content = function(file){
      
      Mutation_subtype = Mutation_subtype_cptac()
      genes = gene_cptac_subtype()
      cancer_type2_subtype = cancer_type2_subtype()
      CPTAC_cohort_subtype = CPTAC_cohort_subtype()
      
      if(input$CPTAC_sur_subtype_res <= 300){
        r = input$CPTAC_sur_subtype_res
      }else{
        r = 300
      }
      if(input$CPTAC_sur_subtype_file == "pdf" & input$CPTAC_sur_subtype_height <= 50){
        h = input$CPTAC_sur_subtype_height
      }else if(input$CPTAC_sur_subtype_file == "pdf" & input$CPTAC_sur_subtype_height > 50){
        h = 50
      }else if(input$CPTAC_sur_subtype_file != "pdf" & input$CPTAC_sur_subtype_height <= 3000){
        h = input$CPTAC_sur_subtype_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sur_subtype_file == "pdf" & input$CPTAC_sur_subtype_width <= 50){
        w = input$CPTAC_sur_subtype_width
      }else if(input$CPTAC_sur_subtype_file == "pdf" & input$CPTAC_sur_subtype_width > 50){
        w = 50
      }else if(input$CPTAC_sur_subtype_file != "pdf" & input$CPTAC_sur_subtype_width <= 3000){
        w = input$CPTAC_sur_subtype_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sur_subtype_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sur_subtype_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sur_subtype_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      tmp_data = CPTAC_survival_subtype_fun_tmpdata_compare(TCGA = CPTAC,cancer_type = cancer_type2_subtype,genes = genes,mut_patient_id = CPTAC_cohort_subtype$mut,wt_patient_id = CPTAC_cohort_subtype$wt)
      print(CPTAC_survival_subtype_fun_compare(tmp_data_compare = tmp_data,genes = genes))
      
      dev.off()
    }
  )
  
  
}