CPTAC_server_pm <- function(input,output,session,CPTAC_pm_para,CPTAC,pathway_database,
                            min.pct_CPTAC_rna_pm,FC_CPTAC_rna_pm,pvalue_CPTAC_rna_pm,min.pct_CPTAC_protein_pm,FC_CPTAC_protein_pm,pvalue_CPTAC_protein_pm,useid){
  
  useidr <- paste(useid,"CPTACpmr",sep = "_")
  useidp <- paste(useid,"CPTACpmp",sep = "_")
  
  CPTAC_cohort_pm = reactive({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    Mut_type_cptac_pm = CPTAC_pm_para()$Mut_type_cptac_pm
    Wild_type_cptac_pm = CPTAC_pm_para()$Wild_type_cptac_pm
    validate(need(any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL),message = FALSE))
    
    CPTAC_cohort_cal_pm(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,Mut_type = Mut_type_cptac_pm,Wild_type = Wild_type_cptac_pm)
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_pm_para()$Mut_type_cptac_pm,CPTAC_pm_para()$Wild_type_cptac_pm)
  
  
  ##################################CPTAC_mutation################################################
  
  output$mut_cptac_pm_uidown1 = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_mut_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_mut_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_mut_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_mut_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_mut_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_mut_pm_file,{
    if(input$CPTAC_mut_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_mut_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_mut_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_mut_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_mut_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_mut_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_mut_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_cptac_pm_uidown2 = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_cm_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_cm_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_cm_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_cm_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_cm_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_cm_pm_file,{
    if(input$CPTAC_cm_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_cm_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_cm_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_cm_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_cm_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_cm_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_cm_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_cptac_pm_uitabledown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      downloadButton('CPTAC_cm_pm_tabdown',label = 'Download Table')
    )
  })
  
  observe({
    gene = CPTAC_pm_para()$gene_cptac_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    
    if(length(CPTAC_cohort_pm$mut) <3 | length(CPTAC_cohort_pm$wt) <3){
      closeAlert(session,"warning_cptac_pm_id")
      createAlert(session, "warning_cptac_pm", "warning_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning_cptac_pm_id")
      createAlert(session, "warning_cptac_pm", "warning_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type2,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_cptac_pm_id")
    }
  })
  
  
  
  Mut_res_cptac_pm =reactive({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = Mutational_Landscape_cptac_pm(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,mut = CPTAC_cohort_pm$mut,wt = CPTAC_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm())
  
  
  
  output$maf4_cptac_pm = renderPlot({
    
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    Mut_res_cptac_pm = Mut_res_cptac_pm()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tl = list(Mcolor)
    names(tl) = gene
    tmp = oncoplot(maf = Mut_res_cptac_pm$tmp_maf,
                   legendFontSize = 2,
                   gene_mar = 8,
                   removeNonMutated = FALSE,
                   legend_height=6,
                   fontSize = 1.2,
                   annotationColor = tl,
                   annotationFontSize = 2,
                   annotationDat = Mut_res_cptac_pm$meta,
                   clinicalFeatures = gene,
                   sortByAnnotation = T,
                   genes = Mut_res_cptac_pm$genes
                   )
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,Mut_res_cptac_pm())
  
  
  output$CPTAC_mut_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_mut_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","mut", '.',"pdf")
      }else if(input$CPTAC_mut_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","mut", '.',"png")
      }else if(input$CPTAC_mut_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","mut", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","mut", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_mut_pm_down")
      gene = CPTAC_pm_para()$gene_cptac_pm
      cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
      Mut_res_cptac_pm = Mut_res_cptac_pm()
      
      if(input$CPTAC_mut_pm_res <= 300){
        r = input$CPTAC_mut_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_mut_pm_file == "pdf" & input$CPTAC_mut_pm_height <= 50){
        h = input$CPTAC_mut_pm_height
      }else if(input$CPTAC_mut_pm_file == "pdf" & input$CPTAC_mut_pm_height > 50){
        h = 50
      }else if(input$CPTAC_mut_pm_file != "pdf" & input$CPTAC_mut_pm_height <= 3000){
        h = input$CPTAC_mut_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_mut_pm_file == "pdf" & input$CPTAC_mut_pm_width <= 50){
        w = input$CPTAC_mut_pm_width
      }else if(input$CPTAC_mut_pm_file == "pdf" & input$CPTAC_mut_pm_width > 50){
        w = 50
      }else if(input$CPTAC_mut_pm_file != "pdf" & input$CPTAC_mut_pm_width <= 3000){
        w = input$CPTAC_mut_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_mut_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_mut_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_mut_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      tl = list(Mcolor)
      names(tl) = gene
      oncoplot(maf = Mut_res_cptac_pm$tmp_maf,
               legendFontSize = 2,
               gene_mar = 8,
               removeNonMutated = FALSE,
               legend_height=6,
               fontSize = 1.2,
               annotationColor = tl,
               annotationFontSize = 2,
               annotationDat = Mut_res_cptac_pm$meta,
               clinicalFeatures = gene,
               sortByAnnotation = T,
               genes = Mut_res_cptac_pm$genes
               )
      
      dev.off()
      shinyjs::enable("CPTAC_mut_pm_down")
    }
  )
  
  output$maf5_cptac_pm = renderPlot({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    Mut_res_cptac_pm = Mut_res_cptac_pm()
    
    num = length(Mut_res_cptac_pm$fvsm$results$Hugo_Symbol)
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
    tmp = oncoplot(maf = Mut_res_cptac_pm$tmp_maf,
                   legendFontSize = 2,
                   gene_mar = 8,
                   removeNonMutated = FALSE,
                   legend_height=6,
                   fontSize = 1.2,
                   annotationColor = tl,
                   annotationFontSize = 2,
                   annotationDat = Mut_res_cptac_pm$meta,
                   clinicalFeatures = gene,
                   sortByAnnotation = T,
                   genes = Mut_res_cptac_pm$fvsm$results$Hugo_Symbol[1:num],
                   keepGeneOrder = T
                   )
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,Mut_res_cptac_pm())
  
  output$CPTAC_cm_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_cm_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","cm", '.',"pdf")
      }else if(input$CPTAC_cm_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","cm", '.',"png")
      }else if(input$CPTAC_cm_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","cm", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","cm", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_cm_pm_down")
      gene = CPTAC_pm_para()$gene_cptac_pm
      cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
      Mut_res_cptac_pm = Mut_res_cptac_pm()
      
      if(input$CPTAC_cm_pm_res <= 300){
        r = input$CPTAC_cm_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_cm_pm_file == "pdf" & input$CPTAC_cm_pm_height <= 50){
        h = input$CPTAC_cm_pm_height
      }else if(input$CPTAC_cm_pm_file == "pdf" & input$CPTAC_cm_pm_height > 50){
        h = 50
      }else if(input$CPTAC_cm_pm_file != "pdf" & input$CPTAC_cm_pm_height <= 3000){
        h = input$CPTAC_cm_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_cm_pm_file == "pdf" & input$CPTAC_cm_pm_width <= 50){
        w = input$CPTAC_cm_pm_width
      }else if(input$CPTAC_cm_pm_file == "pdf" & input$CPTAC_cm_pm_width > 50){
        w = 50
      }else if(input$CPTAC_cm_pm_file != "pdf" & input$CPTAC_cm_pm_width <= 3000){
        w = input$CPTAC_cm_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_cm_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_cm_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_cm_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      num = length(Mut_res_cptac_pm$fvsm$results$Hugo_Symbol)
      if(num > 30){
        num =30
      }
      tl = list(Mcolor)
      names(tl) = gene
      oncoplot(maf = Mut_res_cptac_pm$tmp_maf,
               legendFontSize = 2,
               gene_mar = 8,
               removeNonMutated = FALSE,
               legend_height=6,
               fontSize = 1.2,
               annotationColor = tl,
               annotationFontSize = 2,
               annotationDat = Mut_res_cptac_pm$meta,
               clinicalFeatures = gene,
               sortByAnnotation = T,
               genes = Mut_res_cptac_pm$fvsm$results$Hugo_Symbol[1:num],
               keepGeneOrder = T
               )
      
      dev.off()
      shinyjs::enable("CPTAC_cm_pm_down")
    }
  )
  
  
  output$compare_mutation_cptac_pm = renderReactable({
    Mut_res_cptac_pm = Mut_res_cptac_pm()$fvsm$results
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    Mut_res_cptac_pm$Hugo_Symbol = as.factor(Mut_res_cptac_pm$Hugo_Symbol)
    tmp = reactable( Mut_res_cptac_pm,
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
  }) %>% bindCache(Mut_res_cptac_pm())
  
  output$CPTAC_cm_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      paste0(gene,"_",input$Cancer_type2_pm,".","cm", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_cm_pm_tabdown")
      write.csv(x = Mut_res_cptac_pm()$fvsm$results,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("CPTAC_cm_pm_tabdown")
    }
  )
  
  ###################################################CPTAC immune_infiltration ################################################
  
  output$infr_cptac_pm_uidown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_infr_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_infr_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infr_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infr_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_infr_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_infr_pm_file,{
    if(input$CPTAC_infr_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_infr_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_infr_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_infr_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_infr_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infr_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infr_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$infp_cptac_pm_uidown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_infp_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_infp_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infp_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infp_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_infp_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_infp_pm_file,{
    if(input$CPTAC_infp_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_infp_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_infp_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_infp_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_infp_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infp_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infp_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = CPTAC_pm_para()$gene_cptac_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    
    if(length(CPTAC_cohort_pm$mut) <3 | length(CPTAC_cohort_pm$wt) <3){
      closeAlert(session,"warning2_cptac_pm_id")
      createAlert(session, "warning2_cptac_pm", "warning2_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning2_cptac_pm_id")
      createAlert(session, "warning2_cptac_pm", "warning2_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type2,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning2_cptac_pm_id")
    }
  })
  
  
  
  
  output$immune_infiltration1_CPTAC_rna_pm = renderPlotly({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_infiltration_CPTAC_rna_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
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
    
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm())
  
  
  output$CPTAC_infr_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_infr_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","infr", '.',"pdf")
      }else if(input$CPTAC_infr_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","infr", '.',"png")
      }else if(input$CPTAC_infr_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","infr", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","infr", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_infr_pm_down")
      gene = CPTAC_pm_para()$gene_cptac_pm
      cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
      CPTAC_cohort_pm = CPTAC_cohort_pm()
      
      if(input$CPTAC_infr_pm_res <= 300){
        r = input$CPTAC_infr_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_infr_pm_file == "pdf" & input$CPTAC_infr_pm_height <= 50){
        h = input$CPTAC_infr_pm_height
      }else if(input$CPTAC_infr_pm_file == "pdf" & input$CPTAC_infr_pm_height > 50){
        h = 50
      }else if(input$CPTAC_infr_pm_file != "pdf" & input$CPTAC_infr_pm_height <= 3000){
        h = input$CPTAC_infr_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_infr_pm_file == "pdf" & input$CPTAC_infr_pm_width <= 50){
        w = input$CPTAC_infr_pm_width
      }else if(input$CPTAC_infr_pm_file == "pdf" & input$CPTAC_infr_pm_width > 50){
        w = 50
      }else if(input$CPTAC_infr_pm_file != "pdf" & input$CPTAC_infr_pm_width <= 3000){
        w = input$CPTAC_infr_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_infr_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_infr_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_infr_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(    immune_infiltration_CPTAC_rna_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)     )
      
      dev.off()
      shinyjs::enable("CPTAC_infr_pm_down")
    }
  )
  
  
  output$immune_infiltration1_CPTAC_protein_pm = renderPlotly({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_infiltration_CPTAC_protein_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
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
    
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm())
  
  
  output$CPTAC_infp_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_infp_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","infp", '.',"pdf")
      }else if(input$CPTAC_infp_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","infp", '.',"png")
      }else if(input$CPTAC_infp_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","infp", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","infp", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_infp_pm_down")
      gene = CPTAC_pm_para()$gene_cptac_pm
      cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
      CPTAC_cohort_pm = CPTAC_cohort_pm()
      
      if(input$CPTAC_infp_pm_res <= 300){
        r = input$CPTAC_infp_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_infp_pm_file == "pdf" & input$CPTAC_infp_pm_height <= 50){
        h = input$CPTAC_infp_pm_height
      }else if(input$CPTAC_infp_pm_file == "pdf" & input$CPTAC_infp_pm_height > 50){
        h = 50
      }else if(input$CPTAC_infp_pm_file != "pdf" & input$CPTAC_infp_pm_height <= 3000){
        h = input$CPTAC_infp_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_infp_pm_file == "pdf" & input$CPTAC_infp_pm_width <= 50){
        w = input$CPTAC_infp_pm_width
      }else if(input$CPTAC_infp_pm_file == "pdf" & input$CPTAC_infp_pm_width > 50){
        w = 50
      }else if(input$CPTAC_infp_pm_file != "pdf" & input$CPTAC_infp_pm_width <= 3000){
        w = input$CPTAC_infp_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_infp_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_infp_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_infp_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(    immune_infiltration_CPTAC_protein_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)     )
      
      dev.off()
      shinyjs::enable("CPTAC_infp_pm_down")
    }
  )
  
  
  output$immune_infiltration2_CPTAC_rna_pm = renderPlot({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_infiltration",selection = input$immune_cell_type_CPTAC_pm,outlier = input$outlier2_CPTAC_pm,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 400,height = 500) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm(),input$immune_cell_type_CPTAC_pm,input$outlier2_CPTAC_pm)
  
  output$immune_infiltration2_CPTAC_protein_pm = renderPlot({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_infiltration_protein",selection = input$immune_cell_type_CPTAC_pm,outlier = input$outlier2_CPTAC_pm,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    },width = 400,height = 500) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm(),input$immune_cell_type_CPTAC_pm,input$outlier2_CPTAC_pm)
  
  ###################################################CPTAC immune_signature ################################################
  
  output$sigr_cptac_pm_uidown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_sigr_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sigr_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigr_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigr_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sigr_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sigr_pm_file,{
    if(input$CPTAC_sigr_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sigr_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sigr_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$sigp_cptac_pm_uidown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_sigp_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sigp_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigp_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigp_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sigp_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sigp_pm_file,{
    if(input$CPTAC_sigp_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sigp_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sigp_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observe({
    gene = CPTAC_pm_para()$gene_cptac_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    
    if(length(CPTAC_cohort_pm$mut) <3 | length(CPTAC_cohort_pm$wt) <3){
      closeAlert(session,"warning3_cptac_pm_id")
      createAlert(session, "warning3_cptac_pm", "warning3_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning3_cptac_pm_id")
      createAlert(session, "warning3_cptac_pm", "warning3_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type2,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning3_cptac_pm_id")
    }
  })
  
  
  
  output$immune_signature1_CPTAC_rna_pm = renderPlotly({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_signature_CPTAC_rna_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
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
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm())
  
  
  output$CPTAC_sigr_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_sigr_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","sigr", '.',"pdf")
      }else if(input$CPTAC_sigr_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","sigr", '.',"png")
      }else if(input$CPTAC_sigr_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","sigr", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","sigr", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_sigr_pm_down")
      gene = CPTAC_pm_para()$gene_cptac_pm
      cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
      CPTAC_cohort_pm = CPTAC_cohort_pm()
      
      if(input$CPTAC_sigr_pm_res <= 300){
        r = input$CPTAC_sigr_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_sigr_pm_file == "pdf" & input$CPTAC_sigr_pm_height <= 50){
        h = input$CPTAC_sigr_pm_height
      }else if(input$CPTAC_sigr_pm_file == "pdf" & input$CPTAC_sigr_pm_height > 50){
        h = 50
      }else if(input$CPTAC_sigr_pm_file != "pdf" & input$CPTAC_sigr_pm_height <= 3000){
        h = input$CPTAC_sigr_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sigr_pm_file == "pdf" & input$CPTAC_sigr_pm_width <= 50){
        w = input$CPTAC_sigr_pm_width
      }else if(input$CPTAC_sigr_pm_file == "pdf" & input$CPTAC_sigr_pm_width > 50){
        w = 50
      }else if(input$CPTAC_sigr_pm_file != "pdf" & input$CPTAC_sigr_pm_width <= 3000){
        w = input$CPTAC_sigr_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sigr_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sigr_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sigr_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(    immune_signature_CPTAC_rna_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)     )
      
      dev.off()
      shinyjs::enable("CPTAC_sigr_pm_down")
    }
  )
  
  
  output$immune_signature1_CPTAC_protein_pm = renderPlotly({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_signature_CPTAC_protein_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
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
     
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm())
  
  
  output$CPTAC_sigp_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_sigp_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","sigp", '.',"pdf")
      }else if(input$CPTAC_sigp_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","sigp", '.',"png")
      }else if(input$CPTAC_sigp_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","sigp", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","sigp", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_sigp_pm_down")
      gene = CPTAC_pm_para()$gene_cptac_pm
      cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
      CPTAC_cohort_pm = CPTAC_cohort_pm()
      
      if(input$CPTAC_sigp_pm_res <= 300){
        r = input$CPTAC_sigp_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_sigp_pm_file == "pdf" & input$CPTAC_sigp_pm_height <= 50){
        h = input$CPTAC_sigp_pm_height
      }else if(input$CPTAC_sigp_pm_file == "pdf" & input$CPTAC_sigp_pm_height > 50){
        h = 50
      }else if(input$CPTAC_sigp_pm_file != "pdf" & input$CPTAC_sigp_pm_height <= 3000){
        h = input$CPTAC_sigp_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sigp_pm_file == "pdf" & input$CPTAC_sigp_pm_width <= 50){
        w = input$CPTAC_sigp_pm_width
      }else if(input$CPTAC_sigp_pm_file == "pdf" & input$CPTAC_sigp_pm_width > 50){
        w = 50
      }else if(input$CPTAC_sigp_pm_file != "pdf" & input$CPTAC_sigp_pm_width <= 3000){
        w = input$CPTAC_sigp_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sigp_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sigp_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sigp_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(    immune_signature_CPTAC_protein_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)     )
      
      dev.off()
      shinyjs::enable("CPTAC_sigp_pm_down")
    }
  )
  
  
  output$immune_signature2_CPTAC_rna_pm = renderPlot({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_pathway",selection = input$immune_signature_CPTAC_pm,outlier = input$outlier3_CPTAC_pm,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 400,height = 500) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm(),input$immune_signature_CPTAC_pm,input$outlier3_CPTAC_pm)
  
  output$immune_signature2_CPTAC_protein_pm = renderPlot({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    validate(need(input$immune_signature_CPTAC_pm %in% rownames(CPTAC[[cancer_type2]][["immune_pathway_protein"]]),paste("The pathway",input$immune_signature_CPTAC_pm,"does not exist in",cancer_type2,sep = " ")))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_pathway_protein",selection = input$immune_signature_CPTAC_pm,outlier = input$outlier3_CPTAC_pm,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 400,height = 500) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm(),input$immune_signature_CPTAC_pm,input$outlier3_CPTAC_pm)
  
  ###################################################CPTAC DEG ################################################
  output$DEGr_cptac_pm_uitabledown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_diffr_pm_tabdown',label = 'Download Table')
  })
  
  output$DEGp_cptac_pm_uitabledown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_diffp_pm_tabdown',label = 'Download Table')
  })
  
  observe({
    gene = CPTAC_pm_para()$gene_cptac_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    
    if(length(CPTAC_cohort_pm$mut) <3 | length(CPTAC_cohort_pm$wt) <3){
      closeAlert(session,"warning4_cptac_pm_id")
      createAlert(session, "warning4_cptac_pm", "warning4_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning4_cptac_pm_id")
      createAlert(session, "warning4_cptac_pm", "warning4_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type2,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning4_cptac_pm_id")
    }
  })
  
  
  DEG_table_CPTAC_rna_pm = reactive({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    min.pct_CPTAC_rna_pm = min.pct_CPTAC_rna_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = DEG_CPTAC_rna_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,min.pct = min.pct_CPTAC_rna_pm,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm(),min.pct_CPTAC_rna_pm())
  
  output$DEG_tab_CPTAC_rna_pm = renderReactable({
    DEG_table_CPTAC_rna_pm = DEG_table_CPTAC_rna_pm()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_CPTAC_rna_pm,3),
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
    
    
  }) %>% bindCache(DEG_table_CPTAC_rna_pm())
  
  output$CPTAC_diffr_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      paste0(gene,"_",input$Cancer_type2_pm,".","diffr", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_diffr_pm_tabdown")
      write.csv(x = DEG_table_CPTAC_rna_pm(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("CPTAC_diffr_pm_tabdown")
    }
  )
  
  
  output$volcano_CPTAC_rna_pm = renderPlotly({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    tab = DEG_table_CPTAC_rna_pm()
    FC_CPTAC_rna_pm = FC_CPTAC_rna_pm()
    pvalue_CPTAC_rna_pm = pvalue_CPTAC_rna_pm()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
  
    tmp = volcano_plot(tab = tab,FC = FC_CPTAC_rna_pm,pvalue = pvalue_CPTAC_rna_pm)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,DEG_table_CPTAC_rna_pm(),FC_CPTAC_rna_pm(),pvalue_CPTAC_rna_pm())
  
  #protein
  DEG_table_CPTAC_protein_pm = reactive({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    min.pct_CPTAC_protein_pm = min.pct_CPTAC_protein_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = DEG_CPTAC_protein_pm(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,min.pct = min.pct_CPTAC_protein_pm,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm(),min.pct_CPTAC_protein_pm())
  
  output$DEG_tab_CPTAC_protein_pm = renderReactable({
    DEG_table_CPTAC_protein_pm = DEG_table_CPTAC_protein_pm()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_CPTAC_protein_pm,3),
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
  }) %>% bindCache(DEG_table_CPTAC_protein_pm())
  
  output$CPTAC_diffp_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      paste0(gene,"_",input$Cancer_type2_pm,".","diffp", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_diffp_pm_tabdown")
      write.csv(x = DEG_table_CPTAC_protein_pm(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("CPTAC_diffp_pm_tabdown")
    }
  )
  
  
  output$volcano_CPTAC_protein_pm = renderPlotly({
    gene = CPTAC_pm_para()$gene_cptac_pm
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    tab = DEG_table_CPTAC_protein_pm()
    FC_CPTAC_protein_pm = FC_CPTAC_protein_pm()
    pvalue_CPTAC_protein_pm = pvalue_CPTAC_protein_pm()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_CPTAC_protein_pm,pvalue = pvalue_CPTAC_protein_pm)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,DEG_table_CPTAC_protein_pm(),FC_CPTAC_protein_pm(),pvalue_CPTAC_protein_pm())
  
  
  #################################################CPTAC GSEA######################################################
  
  output$GSEAr_cptac_pm_uitabledown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_gsear_pm_tabdown',label = 'Download Table')
  })
  
  output$GSEAp_cptac_pm_uitabledown = renderUI({
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_gseap_pm_tabdown',label = 'Download Table')
  })
  
  
  observe({
    gene = CPTAC_pm_para()$gene_cptac_pm
    CPTAC_cohort_pm = CPTAC_cohort_pm()
    cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
    
    if(length(CPTAC_cohort_pm$mut) <3 | length(CPTAC_cohort_pm$wt) <3){
      closeAlert(session,"warning5_cptac_pm_id")
      createAlert(session, "warning5_cptac_pm", "warning5_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else if(!any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL)){
      closeAlert(session,"warning5_cptac_pm_id")
      createAlert(session, "warning5_cptac_pm", "warning5_cptac_pm_id", title = "Warning",style = "danger",
                  content = paste("The genes of ",gene,"does not exist in",cancer_type2,sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning5_cptac_pm_id")
    }
  })
  
  status_file_rna <- tempfile()
  write("", status_file_rna)
  status_file_protein <- tempfile()
  write("", status_file_protein)
  
  tmp_gsea_cptac_rna_pm = reactive({
    diff = DEG_table_CPTAC_rna_pm()
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
    
    
    if( length(GSEA_use) < 5 | useidr %in% names(GSEA_use) ){
      
      installr::kill_pid(pid = scan(status_file_rna))
      tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_cptac_rna_pm]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
      write(tmp$get_pid(), status_file_rna)
      
    }else{
      
      closeAlert(session,"warning5_cptac_pm_id")
      createAlert(session, "warning5_cptac_pm", "warning5_cptac_pm_id", title = "Warning",style = "danger",
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
  
  one_loader_rna = reactiveVal(0)
  output$GSEA_tab_cptac_rna_pm = renderReactable({
    if(is.null(tmp_gsea_cptac_rna_pm())){return(NULL)}
    if(tmp_gsea_cptac_rna_pm()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader_rna() == 0){
        
        waiter_show(id = "cptac_pm_rna_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader_rna(one_loader_rna()+1)
        GSEA_use[[useidr]] <<- TRUE
        return(NULL)
        
        
      }else if(one_loader_rna() == 1){
        
        shinyjs::disable(id = "CPTAC_gsear_pm_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "cptac_pm_rna_GSEA_loader")
      
      tmp_gsea_cptac_rna_pm = tmp_gsea_cptac_rna_pm()$get_result()@result
      tmp_gsea_cptac_rna_pm$ID = as.factor(tmp_gsea_cptac_rna_pm$ID)
      tmp_gsea_cptac_rna_pm$Description = as.factor(tmp_gsea_cptac_rna_pm$Description)
      tmp_gsea_cptac_rna_pm$Plot = NA
      tmp_gsea_cptac_rna_pm = tmp_gsea_cptac_rna_pm[,c(ncol(tmp_gsea_cptac_rna_pm),3:(ncol(tmp_gsea_cptac_rna_pm)-2))]
      tmp = reactable( tmp_gsea_cptac_rna_pm,
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
                           Shiny.setInputValue('cptac_pm_rna_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "CPTAC_gsear_pm_tabdown")
      one_loader_rna(0)
      write("", status_file_rna)
      
      GSEA_use[[useidr]] <<- NULL
      return(tmp)
    }
    
    

    
  })
  
  observeEvent(input$cptac_pm_rna_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "CPTAC_gsear_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "CPTAC_gsear_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gsear_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gsear_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "CPTAC_gsear_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_cptac_rna_pm"),type = 1)
    ))
  })
  
  observeEvent(input$CPTAC_gsear_pm_file,{
    if(input$CPTAC_gsear_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_gsear_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_gsear_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$CPTAC_gsear_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      paste0(gene,"_",input$Cancer_type2_pm,".","gsear", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_gsear_pm_tabdown")
      write.csv(x = tmp_gsea_cptac_rna_pm()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("CPTAC_gsear_pm_tabdown")
    }
  )
  
  
  row_num_cptac_rna_pm = eventReactive(input$cptac_pm_rna_GSEA_show,{input$cptac_pm_rna_GSEA_show$index})%>% debounce(500)
  
  output$GSEA_plot_cptac_rna_pm = renderPlot({
    my_gseaplot2(tmp_gsea_cptac_rna_pm()$get_result(),geneSetID = tmp_gsea_cptac_rna_pm()$get_result()@result$ID[row_num_cptac_rna_pm()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_cptac_rna_pm()$get_result()@result$ID[row_num_cptac_rna_pm()])
  })
  
  output$CPTAC_gsear_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_gsear_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","gsear", '.',"pdf")
      }else if(input$CPTAC_gsear_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","gsear", '.',"png")
      }else if(input$CPTAC_gsear_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","gsear", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","gsear", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_gsear_pm_down")
      if(input$CPTAC_gsear_pm_res <= 300){
        r = input$CPTAC_gsear_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_gsear_pm_file == "pdf" & input$CPTAC_gsear_pm_height <= 50){
        h = input$CPTAC_gsear_pm_height
      }else if(input$CPTAC_gsear_pm_file == "pdf" & input$CPTAC_gsear_pm_height > 50){
        h = 50
      }else if(input$CPTAC_gsear_pm_file != "pdf" & input$CPTAC_gsear_pm_height <= 3000){
        h = input$CPTAC_gsear_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_gsear_pm_file == "pdf" & input$CPTAC_gsear_pm_width <= 50){
        w = input$CPTAC_gsear_pm_width
      }else if(input$CPTAC_gsear_pm_file == "pdf" & input$CPTAC_gsear_pm_width > 50){
        w = 50
      }else if(input$CPTAC_gsear_pm_file != "pdf" & input$CPTAC_gsear_pm_width <= 3000){
        w = input$CPTAC_gsear_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_gsear_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_gsear_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_gsear_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(
        my_gseaplot2(tmp_gsea_cptac_rna_pm()$get_result(),geneSetID = tmp_gsea_cptac_rna_pm()$get_result()@result$ID[row_num_cptac_rna_pm()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_cptac_rna_pm()$get_result()@result$ID[row_num_cptac_rna_pm()])
      )
      
      dev.off()
      shinyjs::enable("CPTAC_gsear_pm_down")
    }
  )
  
  
  tmp_gsea_cptac_protein_pm = reactive({
    diff = DEG_table_CPTAC_protein_pm()
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
    
    if( length(GSEA_use) < 5 | useidp %in% names(GSEA_use) ){
      
      installr::kill_pid(pid = scan(status_file_protein))
      tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_cptac_protein_pm]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
      write(tmp$get_pid(), status_file_protein)
      
    }else{
      
      closeAlert(session,"warning5_cptac_pm_id")
      createAlert(session, "warning5_cptac_pm", "warning5_cptac_pm_id", title = "Warning",style = "danger",
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
  
  one_loader_protein = reactiveVal(0)
  output$GSEA_tab_cptac_protein_pm = renderReactable({
    if(is.null(tmp_gsea_cptac_protein_pm())){return(NULL)}
    if(tmp_gsea_cptac_protein_pm()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader_protein() == 0){
        
        waiter_show(id = "cptac_pm_protein_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader_protein(one_loader_protein()+1)
        GSEA_use[[useidp]] <<- TRUE
        return(NULL)
        
        
      }else if(one_loader_protein() == 1){
        
        shinyjs::disable(id = "CPTAC_gseap_pm_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "cptac_pm_protein_GSEA_loader")
      
      tmp_gsea_cptac_protein_pm = tmp_gsea_cptac_protein_pm()$get_result()@result
      tmp_gsea_cptac_protein_pm$ID = as.factor(tmp_gsea_cptac_protein_pm$ID)
      tmp_gsea_cptac_protein_pm$Description = as.factor(tmp_gsea_cptac_protein_pm$Description)
      tmp_gsea_cptac_protein_pm$Plot = NA
      tmp_gsea_cptac_protein_pm = tmp_gsea_cptac_protein_pm[,c(ncol(tmp_gsea_cptac_protein_pm),3:(ncol(tmp_gsea_cptac_protein_pm)-2))]
      tmp = reactable( tmp_gsea_cptac_protein_pm,
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
                           Shiny.setInputValue('cptac_pm_protein_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "CPTAC_gseap_pm_tabdown")
      one_loader_protein(0)
      write("", status_file_protein)
      
      GSEA_use[[useidp]] <<- NULL
      return(tmp)
    }
    
    
    

  })
  
  
  observeEvent(input$cptac_pm_protein_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "CPTAC_gseap_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "CPTAC_gseap_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gseap_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gseap_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "CPTAC_gseap_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_cptac_protein_pm"),type = 1)
    ))
  })
  
  observeEvent(input$CPTAC_gseap_pm_file,{
    if(input$CPTAC_gseap_pm_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_gseap_pm_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_pm_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_pm_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_gseap_pm_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_pm_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_pm_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$CPTAC_gseap_pm_tabdown = downloadHandler(
    
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      paste0(gene,"_",input$Cancer_type2_pm,".","gseap", '.',"csv")
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_gseap_pm_tabdown")
      write.csv(x = tmp_gsea_cptac_protein_pm()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      shinyjs::enable("CPTAC_gseap_pm_tabdown")
    }
  )
  
  
  row_num_cptac_protein_pm = eventReactive(input$cptac_pm_protein_GSEA_show,{input$cptac_pm_protein_GSEA_show$index})%>% debounce(500)
  
  output$GSEA_plot_cptac_protein_pm = renderPlot({
    my_gseaplot2(tmp_gsea_cptac_protein_pm()$get_result(),geneSetID = tmp_gsea_cptac_protein_pm()$get_result()@result$ID[row_num_cptac_protein_pm()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_cptac_protein_pm()$get_result()@result$ID[row_num_cptac_protein_pm()])
  })
  
  
  output$CPTAC_gseap_pm_down = downloadHandler(
    filename = function(){
      gene = CPTAC_pm_para()$gene_cptac_pm
      
      if(input$CPTAC_gseap_pm_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2_pm,".","gseap", '.',"pdf")
      }else if(input$CPTAC_gseap_pm_file == 'png'){
        paste0(gene,"_",input$Cancer_type2_pm,".","gseap", '.',"png")
      }else if(input$CPTAC_gseap_pm_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2_pm,".","gseap", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2_pm,".","gseap", '.',"tiff")
      }
      
    },
    content = function(file){
      shinyjs::disable("CPTAC_gseap_pm_down")
      if(input$CPTAC_gseap_pm_res <= 300){
        r = input$CPTAC_gseap_pm_res
      }else{
        r = 300
      }
      if(input$CPTAC_gseap_pm_file == "pdf" & input$CPTAC_gseap_pm_height <= 50){
        h = input$CPTAC_gseap_pm_height
      }else if(input$CPTAC_gseap_pm_file == "pdf" & input$CPTAC_gseap_pm_height > 50){
        h = 50
      }else if(input$CPTAC_gseap_pm_file != "pdf" & input$CPTAC_gseap_pm_height <= 3000){
        h = input$CPTAC_gseap_pm_height
      }else{
        h = 3000
      }
      if(input$CPTAC_gseap_pm_file == "pdf" & input$CPTAC_gseap_pm_width <= 50){
        w = input$CPTAC_gseap_pm_width
      }else if(input$CPTAC_gseap_pm_file == "pdf" & input$CPTAC_gseap_pm_width > 50){
        w = 50
      }else if(input$CPTAC_gseap_pm_file != "pdf" & input$CPTAC_gseap_pm_width <= 3000){
        w = input$CPTAC_gseap_pm_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_gseap_pm_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_gseap_pm_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_gseap_pm_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      print(
        my_gseaplot2(tmp_gsea_cptac_protein_pm()$get_result(),geneSetID = tmp_gsea_cptac_protein_pm()$get_result()@result$ID[row_num_cptac_protein_pm()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_cptac_protein_pm()$get_result()@result$ID[row_num_cptac_protein_pm()])
      )
      
      dev.off()
      shinyjs::enable("CPTAC_gseap_pm_down")
    }
  )
  
  # ###################CPTAC survival###################
  # 
  # output$sur_cptac_pm_uidown = renderUI({
  #   CPTAC_cohort_pm = CPTAC_cohort_pm()
  #   validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
  #   tagList(
  #     div(selectInput(inputId = "CPTAC_sur_pm_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
  #     div(numericInput(inputId = "CPTAC_sur_pm_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
  #     div(numericInput(inputId = "CPTAC_sur_pm_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
  #     div(numericInput(inputId = "CPTAC_sur_pm_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
  #     div(downloadButton(outputId = "CPTAC_sur_pm_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
  #   )
  # })
  # 
  # observeEvent(input$CPTAC_sur_pm_file,{
  #   if(input$CPTAC_sur_pm_file == "pdf"){
  #     updateNumericInput(session = session,inputId = "CPTAC_sur_pm_res",value = 0,min = 0,max = 0)
  #     updateNumericInput(session = session,inputId = "CPTAC_sur_pm_width",value = 10,min = 1,max = 30,step = 1)
  #     updateNumericInput(session = session,inputId = "CPTAC_sur_pm_height",value = 10,min = 1,max = 30,step = 1)
  #   }else{
  #     updateNumericInput(session = session,inputId = "CPTAC_sur_pm_res",min = 100,max = 300,value = 100,step = 10)
  #     updateNumericInput(session = session,inputId = "CPTAC_sur_pm_width",min = 500,max = 3000,value = 1000,step = 10)
  #     updateNumericInput(session = session,inputId = "CPTAC_sur_pm_height",min = 500,max = 3000,value = 1000,step = 10)
  #   }
  # })
  # 
  # observe({
  #   gene = CPTAC_pm_para()$gene_cptac_pm
  #   CPTAC_cohort_pm = CPTAC_cohort_pm()
  #   cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
  #   
  #   if(length(CPTAC_cohort_pm$mut) <3 | length(CPTAC_cohort_pm$wt) <3){
  #     closeAlert(session,"warning6_cptac_pm_id")
  #     createAlert(session, "warning6_cptac_pm", "warning6_cptac_pm_id", title = "Warning",style = "danger",
  #                 content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
  #   }else if(!any(pathway_list[[gene]] %in% CPTAC[[cancer_type2]][["maf"]]@data$SYMBOL)){
  #     closeAlert(session,"warning6_cptac_pm_id")
  #     createAlert(session, "warning6_cptac_pm", "warning6_cptac_pm_id", title = "Warning",style = "danger",
  #                 content = paste("The genes of ",gene,"does not exist in",cancer_type2,sep = " "), append = FALSE)
  #   }else{
  #     closeAlert(session,"warning6_cptac_pm_id")
  #   }
  # })
  # 
  # 
  # output$CPTAC_survival_pm = renderPlot({
  #   gene = CPTAC_pm_para()$gene_cptac_pm
  #   cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
  #   CPTAC_cohort_pm = CPTAC_cohort_pm()
  #   validate(need(length(CPTAC_cohort_pm$mut) >=3 & length(CPTAC_cohort_pm$wt) >= 3,message = FALSE))
  #   
  #   shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
  #   shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
  #   shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
  #   shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
  #   shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
  #   shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
  #   shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
  #   
  #   
  #   tmp = CPTAC_survival_pm(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)
  #   
  #   shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
  #   shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
  #   shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
  #   shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
  #   shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
  #   shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
  #   shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
  #   return(tmp)
  #   
  # },width = 600,height = 600) %>% bindCache(CPTAC_pm_para()$cancer_type2_pm,CPTAC_pm_para()$gene_cptac_pm,CPTAC_cohort_pm())
  # 
  # output$CPTAC_sur_pm_down = downloadHandler(
  #   filename = function(){
  #     gene = CPTAC_pm_para()$gene_cptac_pm
  #     
  #     if(input$CPTAC_sur_pm_file == 'pdf'){
  #       paste0(gene,"_",input$Cancer_type2_pm,".","sur", '.',"pdf")
  #     }else if(input$CPTAC_sur_pm_file == 'png'){
  #       paste0(gene,"_",input$Cancer_type2_pm,".","sur", '.',"png")
  #     }else if(input$CPTAC_sur_pm_file == 'jpeg'){
  #       paste0(gene,"_",input$Cancer_type2_pm,".","sur", '.',"jpeg")
  #     }else{
  #       paste0(gene,"_",input$Cancer_type2_pm,".","sur", '.',"tiff")
  #     }
  #     
  #   },
  #   content = function(file){
  #     shinyjs::disable("CPTAC_sur_pm_down")
  #     gene = CPTAC_pm_para()$gene_cptac_pm
  #     cancer_type2 = CPTAC_pm_para()$cancer_type2_pm
  #     CPTAC_cohort_pm = CPTAC_cohort_pm()
  #     
  #     if(input$CPTAC_sur_pm_res <= 300){
  #       r = input$CPTAC_sur_pm_res
  #     }else{
  #       r = 300
  #     }
  #     if(input$CPTAC_sur_pm_file == "pdf" & input$CPTAC_sur_pm_height <= 50){
  #       h = input$CPTAC_sur_pm_height
  #     }else if(input$CPTAC_sur_pm_file == "pdf" & input$CPTAC_sur_pm_height > 50){
  #       h = 50
  #     }else if(input$CPTAC_sur_pm_file != "pdf" & input$CPTAC_sur_pm_height <= 3000){
  #       h = input$CPTAC_sur_pm_height
  #     }else{
  #       h = 3000
  #     }
  #     if(input$CPTAC_sur_pm_file == "pdf" & input$CPTAC_sur_pm_width <= 50){
  #       w = input$CPTAC_sur_pm_width
  #     }else if(input$CPTAC_sur_pm_file == "pdf" & input$CPTAC_sur_pm_width > 50){
  #       w = 50
  #     }else if(input$CPTAC_sur_pm_file != "pdf" & input$CPTAC_sur_pm_width <= 3000){
  #       w = input$CPTAC_sur_pm_width
  #     }else{
  #       w = 3000
  #     }
  #     
  #     if(input$CPTAC_sur_pm_file == 'pdf'){
  #       pdf(file = file,width = w,height = h)
  #     }else if(input$CPTAC_sur_pm_file == 'png'){
  #       png(file = file,width = w,height = h,res = r)
  #     }else if(input$CPTAC_sur_pm_file == 'jpeg'){
  #       jpeg(file = file,width = w,height = h,res = r)
  #     }else{
  #       tiff(file = file,width = w,height = h,res = r)
  #     }
  # 
  #     
  #     print(    CPTAC_survival_pm(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,mut_patient_id = CPTAC_cohort_pm$mut,wt_patient_id = CPTAC_cohort_pm$wt)     )
  #     
  #     dev.off()
  #     shinyjs::enable("CPTAC_sur_pm_down")
  #   }
  # )
  # 
  # 
}