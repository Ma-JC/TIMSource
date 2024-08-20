TCGA_server_pm <- function(input,output,session,TCGA_pm_para,TCGA,pathway_database,min.pct_pm,FC_pm,pvalue_pm,useid){
  
  useid <- paste(useid,"TCGApm",sep = "_")
  
  TCGA_cohort_pm = reactive({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    Mut_type_tcga_pm = TCGA_pm_para()$Mut_type_tcga_pm
    Wild_type_tcga_pm = TCGA_pm_para()$Wild_type_tcga_pm
    # validate(need(any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL),paste("The genes of ",gene,"does not exist in",cancer_type,sep = " ")))
    validate(need(any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL),message = FALSE))
    TCGA_cohort_cal_pm(TCGA = TCGA,cancer_type = cancer_type,gene = gene,Mut_type = Mut_type_tcga_pm,Wild_type = Wild_type_tcga_pm)
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_pm_para()$Mut_type_tcga_pm,TCGA_pm_para()$Wild_type_tcga_pm)
  
  ##################################TCGA_mutation################################################
  
  output$mut_tcga_pm_uidown1 = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_mut_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_mut_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_mut_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_mut_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_mut_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_mut_pm_file,{
    if(input$TCGA_mut_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_mut_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_mut_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_mut_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_mut_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_mut_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_mut_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_tcga_pm_uidown2 = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_cm_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_cm_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_cm_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_cm_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_cm_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_cm_pm_file,{
    if(input$TCGA_cm_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_cm_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_cm_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_cm_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_cm_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_cm_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_cm_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_cptac_pm_uitabledown = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      downloadButton('TCGA_cm_pm_tabdown',label = 'Download Table')
    )
  })
  
  observe({
    gene = TCGA_pm_para()$gene_tcga_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    cancer_type = TCGA_pm_para()$cancer_type_pm
    
    if(length(TCGA_cohort_pm$mut) <3 | length(TCGA_cohort_pm$wt) <3){
      closeAlert(session,"warning_tcga_pm_id")
      createAlert(session, "warning_tcga_pm", "warning_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning_tcga_pm_id")
      createAlert(session, "warning_tcga_pm", "warning_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_tcga_pm_id")
    }
  })
  
  Mut_res_pm =reactive({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = Mutational_Landscape_pm(TCGA = TCGA,cancer_type = cancer_type,gene = gene,mut = TCGA_cohort_pm$mut,wt = TCGA_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm())
  
  
  
  output$maf4_pm = renderPlot({

    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    Mut_res_pm = Mut_res_pm()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tl = list(Mcolor)
    names(tl) = gene
    tmp = oncoplot(maf = Mut_res_pm$tmp_maf,
                   legendFontSize = 2,
                   gene_mar = 8,
                   removeNonMutated = FALSE,
                   legend_height=6,
                   fontSize = 1.2,
                   annotationColor = tl,
                   annotationFontSize = 2,
                   annotationDat = Mut_res_pm$meta,
                   clinicalFeatures = gene,
                   sortByAnnotation = T,
                   genes = Mut_res_pm$genes
                   )
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,Mut_res_pm())
  
  output$TCGA_mut_pm_down = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      
      if(input$TCGA_mut_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type_pm,".","mut", '.',"pdf")
      }else if(input$TCGA_mut_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type_pm,".","mut", '.',"png")
      }else if(input$TCGA_mut_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type_pm,".","mut", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type_pm,".","mut", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("TCGA_mut_pm_down")
      gene = TCGA_pm_para()$gene_tcga_pm
      cancer_type = TCGA_pm_para()$cancer_type_pm
      Mut_res_pm = Mut_res_pm()
      
      if(input$TCGA_mut_pm_res <= 300){
        r = input$TCGA_mut_pm_res
      }else{
        r = 300
      }
      if(input$TCGA_mut_pm_file == "pdf" & input$TCGA_mut_pm_height <= 50){
        h = input$TCGA_mut_pm_height
      }else if(input$TCGA_mut_pm_file == "pdf" & input$TCGA_mut_pm_height > 50){
        h = 50
      }else if(input$TCGA_mut_pm_file != "pdf" & input$TCGA_mut_pm_height <= 3000){
        h = input$TCGA_mut_pm_height
      }else{
        h = 3000
      }
      if(input$TCGA_mut_pm_file == "pdf" & input$TCGA_mut_pm_width <= 50){
        w = input$TCGA_mut_pm_width
      }else if(input$TCGA_mut_pm_file == "pdf" & input$TCGA_mut_pm_width > 50){
        w = 50
      }else if(input$TCGA_mut_pm_file != "pdf" & input$TCGA_mut_pm_width <= 3000){
        w = input$TCGA_mut_pm_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_mut_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_mut_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_mut_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tl = list(Mcolor)
      names(tl) = gene
      oncoplot(maf = Mut_res_pm$tmp_maf,
               legendFontSize = 2,
               gene_mar = 8,
               removeNonMutated = FALSE,
               legend_height=6,
               fontSize = 1.2,
               annotationColor = tl,
               annotationFontSize = 2,
               annotationDat = Mut_res_pm$meta,
               clinicalFeatures = gene,
               sortByAnnotation = T,
               genes = Mut_res_pm$genes)

      dev.off()
      shinyjs::enable("TCGA_mut_pm_down")
    }
  )

  output$maf5_pm = renderPlot({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    Mut_res_pm = Mut_res_pm()
    num = length(Mut_res_pm$fvsm$results$Hugo_Symbol)
    if(num > 30){
      num =30
    }
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tl = list(Mcolor)
    names(tl) = gene
    tmp = oncoplot(maf = Mut_res_pm$tmp_maf,
                   legendFontSize = 2,
                   gene_mar = 8,
                   removeNonMutated = FALSE,
                   legend_height=6,
                   fontSize = 1.2,
                   annotationColor = tl,
                   annotationFontSize = 2,
                   annotationDat = Mut_res_pm$meta,
                   clinicalFeatures = gene,
                   sortByAnnotation = T,
                   genes = Mut_res_pm$fvsm$results$Hugo_Symbol[1:num],
                   keepGeneOrder = T)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,Mut_res_pm())
  
  output$TCGA_cm_pm_down = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      
      if(input$TCGA_cm_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type_pm,".","cm", '.',"pdf")
      }else if(input$TCGA_cm_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type_pm,".","cm", '.',"png")
      }else if(input$TCGA_cm_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type_pm,".","cm", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type_pm,".","cm", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("TCGA_cm_pm_down")
      gene = TCGA_pm_para()$gene_tcga_pm
      cancer_type = TCGA_pm_para()$cancer_type_pm
      Mut_res_pm = Mut_res_pm()
      
      if(input$TCGA_cm_pm_res <= 300){
        r = input$TCGA_cm_pm_res
      }else{
        r = 300
      }
      if(input$TCGA_cm_pm_file == "pdf" & input$TCGA_cm_pm_height <= 50){
        h = input$TCGA_cm_pm_height
      }else if(input$TCGA_cm_pm_file == "pdf" & input$TCGA_cm_pm_height > 50){
        h = 50
      }else if(input$TCGA_cm_pm_file != "pdf" & input$TCGA_cm_pm_height <= 3000){
        h = input$TCGA_cm_pm_height
      }else{
        h = 3000
      }
      if(input$TCGA_cm_pm_file == "pdf" & input$TCGA_cm_pm_width <= 50){
        w = input$TCGA_cm_pm_width
      }else if(input$TCGA_cm_pm_file == "pdf" & input$TCGA_cm_pm_width > 50){
        w = 50
      }else if(input$TCGA_cm_pm_file != "pdf" & input$TCGA_cm_pm_width <= 3000){
        w = input$TCGA_cm_pm_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_cm_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_cm_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_cm_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      num = length(Mut_res_pm$fvsm$results$Hugo_Symbol)
      if(num > 30){
        num =30
      }
      tl = list(Mcolor)
      names(tl) = gene
      oncoplot(maf = Mut_res_pm$tmp_maf,
               legendFontSize = 2,
               gene_mar = 8,
               removeNonMutated = FALSE,
               legend_height=6,
               fontSize = 1.2,
               annotationColor = tl,
               annotationFontSize = 2,
               annotationDat = Mut_res_pm$meta,
               clinicalFeatures = gene,
               sortByAnnotation = T,
               genes = Mut_res_pm$fvsm$results$Hugo_Symbol[1:num],
               legendFontSize = 2,
               keepGeneOrder = T)
      
      dev.off()
      shinyjs::enable("TCGA_cm_pm_down")
    }
  )
  
  output$compare_mutation_pm = renderReactable({
    Mut_res_pm = Mut_res_pm()$fvsm$results
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    Mut_res_pm$Hugo_Symbol = as.factor(Mut_res_pm$Hugo_Symbol)
    tmp = reactable( Mut_res_pm,
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
    
  }) %>% bindCache(Mut_res_pm())
  
  output$TCGA_cm_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      paste0(gene,"_",input$Cancer_type_pm,".","cm", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("TCGA_cm_pm_tabdown")
      write.csv(x = Mut_res_pm()$fvsm$results,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("TCGA_cm_pm_tabdown")
    }
  )
  ###################################################TCGA immune_infiltration ################################################
  
  output$inf_tcga_pm_uidown = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_inf_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_inf_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_inf_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_inf_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_inf_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_inf_pm_file,{
    if(input$TCGA_inf_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_inf_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_inf_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_inf_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_inf_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_inf_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_inf_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = TCGA_pm_para()$gene_tcga_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    cancer_type = TCGA_pm_para()$cancer_type_pm
    
    if(length(TCGA_cohort_pm$mut) <3 | length(TCGA_cohort_pm$wt) <3){
      closeAlert(session,"warning2_tcga_pm_id")
      createAlert(session, "warning2_tcga_pm", "warning2_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning2_tcga_pm_id")
      createAlert(session, "warning2_tcga_pm", "warning2_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning2_tcga_pm_id")
    }
  })
  
  
  output$immune_infiltration1_pm = renderPlotly({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_infiltration_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)
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
    
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm())
  
  output$TCGA_inf_pm_down = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      
      if(input$TCGA_inf_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type_pm,".","inf", '.',"pdf")
      }else if(input$TCGA_inf_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type_pm,".","inf", '.',"png")
      }else if(input$TCGA_inf_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type_pm,".","inf", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type_pm,".","inf", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("TCGA_inf_pm_down")
      gene = TCGA_pm_para()$gene_tcga_pm
      cancer_type = TCGA_pm_para()$cancer_type_pm
      TCGA_cohort_pm = TCGA_cohort_pm()
      
      if(input$TCGA_inf_pm_res <= 300){
        r = input$TCGA_inf_pm_res
      }else{
        r = 300
      }
      if(input$TCGA_inf_pm_file == "pdf" & input$TCGA_inf_pm_height <= 50){
        h = input$TCGA_inf_pm_height
      }else if(input$TCGA_inf_pm_file == "pdf" & input$TCGA_inf_pm_height > 50){
        h = 50
      }else if(input$TCGA_inf_pm_file != "pdf" & input$TCGA_inf_pm_height <= 3000){
        h = input$TCGA_inf_pm_height
      }else{
        h = 3000
      }
      if(input$TCGA_inf_pm_file == "pdf" & input$TCGA_inf_pm_width <= 50){
        w = input$TCGA_inf_pm_width
      }else if(input$TCGA_inf_pm_file == "pdf" & input$TCGA_inf_pm_width > 50){
        w = 50
      }else if(input$TCGA_inf_pm_file != "pdf" & input$TCGA_inf_pm_width <= 3000){
        w = input$TCGA_inf_pm_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_inf_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_inf_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_inf_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(    immune_infiltration_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)  )
      
      dev.off()
      shinyjs::enable("TCGA_inf_pm_down")
    }
  )
  
  output$immune_infiltration2_pm = renderPlot({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,immune_module = "immune_infiltration",selection = input$immune_cell_type_pm,outlier = input$outlier2_pm,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 600,height = 600) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm(),input$immune_cell_type_pm,input$outlier2_pm)
  ###################################################TCGA immune_signature ################################################
  
  output$sig_tcga_pm_uidown = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_sig_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_sig_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sig_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sig_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_sig_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_sig_pm_file,{
    if(input$TCGA_sig_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_sig_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_sig_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_sig_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_sig_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sig_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sig_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = TCGA_pm_para()$gene_tcga_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    cancer_type = TCGA_pm_para()$cancer_type_pm
    
    if(length(TCGA_cohort_pm$mut) <3 | length(TCGA_cohort_pm$wt) <3){
      closeAlert(session,"warning3_tcga_pm_id")
      createAlert(session, "warning3_tcga_pm", "warning3_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning3_tcga_pm_id")
      createAlert(session, "warning3_tcga_pm", "warning3_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning3_tcga_pm_id")
    }
  })
  
  output$immune_signature1_pm = renderPlotly({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_signature_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)
    tmp = ggplotly(p) %>% 
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0), height = 800) %>% 
            rangeslider(start = 0,end = 9) %>% plotly::config(displayModeBar = FALSE)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm())
  
  
  output$TCGA_sig_pm_down = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      
      if(input$TCGA_sig_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type_pm,".","sig", '.',"pdf")
      }else if(input$TCGA_sig_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type_pm,".","sig", '.',"png")
      }else if(input$TCGA_sig_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type_pm,".","sig", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type_pm,".","sig", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("TCGA_sig_pm_down")
      gene = TCGA_pm_para()$gene_tcga_pm
      cancer_type = TCGA_pm_para()$cancer_type_pm
      TCGA_cohort_pm = TCGA_cohort_pm()
      
      
      if(input$TCGA_sig_pm_res <= 300){
        r = input$TCGA_sig_pm_res
      }else{
        r = 300
      }
      if(input$TCGA_sig_pm_file == "pdf" & input$TCGA_sig_pm_height <= 50){
        h = input$TCGA_sig_pm_height
      }else if(input$TCGA_sig_pm_file == "pdf" & input$TCGA_sig_pm_height > 50){
        h = 50
      }else if(input$TCGA_sig_pm_file != "pdf" & input$TCGA_sig_pm_height <= 3000){
        h = input$TCGA_sig_pm_height
      }else{
        h = 3000
      }
      if(input$TCGA_sig_pm_file == "pdf" & input$TCGA_sig_pm_width <= 50){
        w = input$TCGA_sig_pm_width
      }else if(input$TCGA_sig_pm_file == "pdf" & input$TCGA_sig_pm_width > 50){
        w = 50
      }else if(input$TCGA_sig_pm_file != "pdf" & input$TCGA_sig_pm_width <= 3000){
        w = input$TCGA_sig_pm_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_sig_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_sig_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_sig_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(immune_signature_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)  )
      
      dev.off()
      shinyjs::enable("TCGA_sig_pm_down")
    }
  )
  
  output$immune_signature2_pm = renderPlot({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,immune_module = "immune_pathway",selection = input$immune_signature_pm,outlier = input$outlier3_pm,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 600,height = 600) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm(),input$immune_signature_pm,input$outlier3_pm)
  
  
  ###################################################TCGA DEG ################################################
  
  output$DEG_tcga_pm_uitabledown = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    downloadButton('TCGA_diff_pm_tabdown',label = 'Download Table')
  })
  
  observe({
    gene = TCGA_pm_para()$gene_tcga_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    cancer_type = TCGA_pm_para()$cancer_type_pm
    
    if(length(TCGA_cohort_pm$mut) <3 | length(TCGA_cohort_pm$wt) <3){
      closeAlert(session,"warning4_tcga_pm_id")
      createAlert(session, "warning4_tcga_pm", "warning4_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning4_tcga_pm_id")
      createAlert(session, "warning4_tcga_pm", "warning4_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning4_tcga_pm_id")
    }
  })
  
  DEG_table_pm = reactive({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    min.pct_pm = min.pct_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = DEG_pm(TCGA = TCGA,gene = gene,cancer_type = cancer_type,min.pct = min.pct_pm,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm(),min.pct_pm())
  
  output$DEG_tab_pm = renderReactable({
    DEG_table_pm = DEG_table_pm()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_pm,3),
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
  }) %>% bindCache(DEG_table_pm())
  
  
  output$TCGA_diff_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      paste0(gene,"_",input$Cancer_type_pm,".","diff", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("TCGA_diff_pm_tabdown")
      write.csv(x = DEG_table_pm(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("TCGA_diff_pm_tabdown")
    }
  )
  
  
  output$volcano_pm = renderPlotly({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    tab = DEG_table_pm()
    FC_pm = FC_pm()
    pvalue_pm = pvalue_pm()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_pm,pvalue = pvalue_pm)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,DEG_table_pm(),FC_pm(),pvalue_pm())
  
  
  
  ###################################################TCGA GSEA################################################
  
  output$GSEA_tcga_pm_uitabledown = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    downloadButton('TCGA_gsea_pm_tabdown',label = 'Download Table')
  })
  
  observe({
    gene = TCGA_pm_para()$gene_tcga_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    cancer_type = TCGA_pm_para()$cancer_type_pm
    
    if(length(TCGA_cohort_pm$mut) <3 | length(TCGA_cohort_pm$wt) <3){
      closeAlert(session,"warning5_tcga_pm_id")
      createAlert(session, "warning5_tcga_pm", "warning5_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning5_tcga_pm_id")
      createAlert(session, "warning5_tcga_pm", "warning5_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning5_tcga_pm_id")
    }
  })
  
  status_file <- tempfile()
  write("", status_file)
  
  tmp_gsea_tcga_pm = reactive({
    diff = DEG_table_pm()
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
      tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_tcga_pm]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
      write(tmp$get_pid(), status_file)
      
    }else{
      
      closeAlert(session,"warning5_tcga_pm_id")
      createAlert(session, "warning5_tcga_pm", "warning5_tcga_pm_id", title = "Warning",style = "danger",
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
  output$GSEA_tab_tcga_pm = renderReactable({
    if(is.null(tmp_gsea_tcga_pm())){return(NULL)}
    if(tmp_gsea_tcga_pm()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader() == 0){
        
        waiter_show(id = "tcga_pm_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader(one_loader()+1)
        GSEA_use[[useid]] <<- TRUE
        return(NULL)
        
        
      }else if(one_loader() == 1){
        
        shinyjs::disable(id = "TCGA_gsea_pm_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "tcga_pm_GSEA_loader")
      
      tmp_gsea_tcga_pm = tmp_gsea_tcga_pm()$get_result()@result
      tmp_gsea_tcga_pm$ID = as.factor(tmp_gsea_tcga_pm$ID)
      tmp_gsea_tcga_pm$Description = as.factor(tmp_gsea_tcga_pm$Description)
      tmp_gsea_tcga_pm$Plot = NA
      tmp_gsea_tcga_pm = tmp_gsea_tcga_pm[,c(ncol(tmp_gsea_tcga_pm),3:(ncol(tmp_gsea_tcga_pm)-2))]
      tmp = reactable( tmp_gsea_tcga_pm,
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
                           Shiny.setInputValue('tcga_pm_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      shinyjs::enable(id = "TCGA_gsea_pm_tabdown")
      one_loader(0)
      write("", status_file)
      
      GSEA_use[[useid]] <<- NULL
      return(tmp)
    }
    
    
     return(tmp)
   
  })
  
  observeEvent(input$tcga_pm_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "TCGA_gsea_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "TCGA_gsea_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "TCGA_gsea_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "TCGA_gsea_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "TCGA_gsea_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_tcga_pm"),type = 1)
    ))
  })
  
  observeEvent(input$TCGA_gsea_pm_file,{
    if(input$TCGA_gsea_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_gsea_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_gsea_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_gsea_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_gsea_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_gsea_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_gsea_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$TCGA_gsea_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      paste0(gene,"_",input$Cancer_type_pm,".","gsea", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("TCGA_gsea_pm_tabdown")
      write.csv(x = tmp_gsea_tcga_pm()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("TCGA_gsea_pm_tabdown")
    }
  )
  
  
  row_num_tcga_pm = eventReactive(input$tcga_pm_GSEA_show,{input$tcga_pm_GSEA_show$index})%>% debounce(500)
  
  output$GSEA_plot_tcga_pm = renderPlot({
    my_gseaplot2(tmp_gsea_tcga_pm()$get_result(),geneSetID = tmp_gsea_tcga_pm()$get_result()@result$ID[row_num_tcga_pm()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_tcga_pm()$get_result()@result$ID[row_num_tcga_pm()])
  })
  
  output$TCGA_gsea_pm_down = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      
      if(input$TCGA_gsea_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type_pm,".","gsea", '.',"pdf")
      }else if(input$TCGA_gsea_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type_pm,".","gsea", '.',"png")
      }else if(input$TCGA_gsea_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type_pm,".","gsea", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type_pm,".","gsea", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("TCGA_gsea_pm_down")
      if(input$TCGA_gsea_pm_res <= 300){
        r = input$TCGA_gsea_pm_res
      }else{
        r = 300
      }
      if(input$TCGA_gsea_pm_file == "pdf" & input$TCGA_gsea_pm_height <= 50){
        h = input$TCGA_gsea_pm_height
      }else if(input$TCGA_gsea_pm_file == "pdf" & input$TCGA_gsea_pm_height > 50){
        h = 50
      }else if(input$TCGA_gsea_pm_file != "pdf" & input$TCGA_gsea_pm_height <= 3000){
        h = input$TCGA_gsea_pm_height
      }else{
        h = 3000
      }
      if(input$TCGA_gsea_pm_file == "pdf" & input$TCGA_gsea_pm_width <= 50){
        w = input$TCGA_gsea_pm_width
      }else if(input$TCGA_gsea_pm_file == "pdf" & input$TCGA_gsea_pm_width > 50){
        w = 50
      }else if(input$TCGA_gsea_pm_file != "pdf" & input$TCGA_gsea_pm_width <= 3000){
        w = input$TCGA_gsea_pm_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_gsea_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_gsea_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_gsea_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(
        my_gseaplot2(tmp_gsea_tcga_pm()$get_result(),geneSetID = tmp_gsea_tcga_pm()$get_result()@result$ID[row_num_tcga_pm()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_tcga_pm()$get_result()@result$ID[row_num_tcga_pm()])
      )
      
      dev.off()
      shinyjs::enable("TCGA_gsea_pm_down")
    }
  )
  
  
  ###################TCGA survival###################
  
  
  output$sur_tcga_pm_uidown = renderUI({
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_sur_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_sur_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sur_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sur_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_sur_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_sur_pm_file,{
    if(input$TCGA_sur_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_sur_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_sur_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_sur_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_sur_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sur_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sur_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = TCGA_pm_para()$gene_tcga_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    cancer_type = TCGA_pm_para()$cancer_type_pm
    
    if(length(TCGA_cohort_pm$mut) <3 | length(TCGA_cohort_pm$wt) <3){
      closeAlert(session,"warning6_tcga_pm_id")
      createAlert(session, "warning6_tcga_pm", "warning6_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% TCGA[[cancer_type]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning6_tcga_pm_id")
      createAlert(session, "warning6_tcga_pm", "warning6_tcga_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning6_tcga_pm_id")
    }
  })
  
  width = reactive({    if(TCGA_pm_para()$cancer_type_pm == "LAML"){500}else{1600}})
  output$TCGA_survival_pm = renderPlot({
    gene = TCGA_pm_para()$gene_tcga_pm
    cancer_type = TCGA_pm_para()$cancer_type_pm
    TCGA_cohort_pm = TCGA_cohort_pm()
    validate(need(length(TCGA_cohort_pm$mut) >=3 & length(TCGA_cohort_pm$wt) >= 3,message = FALSE))
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = TCGA_survival_pm(TCGA = TCGA,cancer_type = cancer_type,gene = gene,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    
  },width = width,height = 600)  %>% bindCache(TCGA_pm_para()$cancer_type_pm,TCGA_pm_para()$gene_tcga_pm,TCGA_cohort_pm())
  
  
  output$TCGA_sur_pm_down = downloadHandler(
    
    filename = function(){
      gene = TCGA_pm_para()$gene_tcga_pm
      
      if(input$TCGA_sur_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type_pm,".","sur", '.',"pdf")
      }else if(input$TCGA_sur_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type_pm,".","sur", '.',"png")
      }else if(input$TCGA_sur_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type_pm,".","sur", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type_pm,".","sur", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("TCGA_sur_pm_down")
      gene = TCGA_pm_para()$gene_tcga_pm
      cancer_type = TCGA_pm_para()$cancer_type_pm
      TCGA_cohort_pm = TCGA_cohort_pm()
      
      if(input$TCGA_sur_pm_res <= 300){
        r = input$TCGA_sur_pm_res
      }else{
        r = 300
      }
      if(input$TCGA_sur_pm_file == "pdf" & input$TCGA_sur_pm_height <= 50){
        h = input$TCGA_sur_pm_height
      }else if(input$TCGA_sur_pm_file == "pdf" & input$TCGA_sur_pm_height > 50){
        h = 50
      }else if(input$TCGA_sur_pm_file != "pdf" & input$TCGA_sur_pm_height <= 3000){
        h = input$TCGA_sur_pm_height
      }else{
        h = 3000
      }
      if(input$TCGA_sur_pm_file == "pdf" & input$TCGA_sur_pm_width <= 50){
        w = input$TCGA_sur_pm_width
      }else if(input$TCGA_sur_pm_file == "pdf" & input$TCGA_sur_pm_width > 50){
        w = 50
      }else if(input$TCGA_sur_pm_file != "pdf" & input$TCGA_sur_pm_width <= 3000){
        w = input$TCGA_sur_pm_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_sur_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_sur_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_sur_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      print(  TCGA_survival_pm(TCGA = TCGA,cancer_type = cancer_type,gene = gene,mut_patient_id = TCGA_cohort_pm$mut,wt_patient_id = TCGA_cohort_pm$wt)   )
      
      dev.off()
      shinyjs::enable("TCGA_sur_pm_down")
    }
  )
  
  
  
}

