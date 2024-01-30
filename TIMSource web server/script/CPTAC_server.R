CPTAC_server <- function(input,output,session,gene_cptac,cancer_type2,CPTAC,pathway_database,Mut_type_cptac_single,Wild_type_cptac_single,
                         min.pct_CPTAC_rna,FC_CPTAC_rna,pvalue_CPTAC_rna,min.pct_CPTAC_protein,FC_CPTAC_protein,pvalue_CPTAC_protein){
  
  
  CPTAC_cohort = reactive({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    Mut_type_cptac_single = Mut_type_cptac_single()
    Wild_type_cptac_single = Wild_type_cptac_single()
    
    CPTAC_cohort_cal(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,Mut_type = Mut_type_cptac_single,Wild_type = Wild_type_cptac_single)
  }) %>% bindCache(cancer_type2(),gene_cptac(),Mut_type_cptac_single(),Wild_type_cptac_single())
  
  
  
  ###################################CPTAC_overview################################################
  
  output$over_cptac_single_uidown1 = renderUI({
    cancer_type2 = cancer_type2()
    tagList(
            div(selectInput(inputId = "CPTAC_plot1_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
            div(numericInput(inputId = "CPTAC_plot1_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(numericInput(inputId = "CPTAC_plot1_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(numericInput(inputId = "CPTAC_plot1_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(downloadButton(outputId = "CPTAC_plot1_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;"))
  })
  
  observeEvent(input$CPTAC_plot1_single_file,{
    if(input$CPTAC_plot1_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_plot1_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_plot1_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_plot1_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_plot1_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_plot1_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_plot1_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$over_cptac_single_uidown2 = renderUI({
    cancer_type2 = cancer_type2()
    tagList(
            div(selectInput(inputId = "CPTAC_plot2_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
            div(numericInput(inputId = "CPTAC_plot2_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(numericInput(inputId = "CPTAC_plot2_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(numericInput(inputId = "CPTAC_plot2_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(downloadButton(outputId = "CPTAC_plot2_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;"))
  })
  
  observeEvent(input$CPTAC_plot2_single_file,{
    if(input$CPTAC_plot2_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_plot2_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_plot2_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_plot2_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_plot2_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_plot2_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_plot2_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$over_cptac_single_uidown3 = renderUI({
    cancer_type2 = cancer_type2()
    tagList(
            div(selectInput(inputId = "CPTAC_plot3_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
            div(numericInput(inputId = "CPTAC_plot3_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(numericInput(inputId = "CPTAC_plot3_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(numericInput(inputId = "CPTAC_plot3_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
            div(downloadButton(outputId = "CPTAC_plot3_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;"))
  })
  
  observeEvent(input$CPTAC_plot3_single_file,{
    if(input$CPTAC_plot3_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_plot3_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_plot3_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_plot3_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_plot3_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_plot3_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_plot3_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  

  output$maf1_CPTAC = renderPlot({
    cancer_type2 = cancer_type2()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = plotmafSummary(maf = CPTAC[[cancer_type2]][["maf"]], rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    
  }) %>% bindCache(cancer_type2())
  
  output$CPTAC_plot1_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_plot1_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","plot1", '.',"pdf")
      }else if(input$CPTAC_plot1_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","plot1", '.',"png")
      }else if(input$CPTAC_plot1_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","plot1", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","plot1", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type2 = cancer_type2()
      
      if(input$CPTAC_plot1_single_res <= 300){
        r = input$CPTAC_plot1_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_plot1_single_file == "pdf" & input$CPTAC_plot1_single_height <= 50){
        h = input$CPTAC_plot1_single_height
      }else if(input$CPTAC_plot1_single_file == "pdf" & input$CPTAC_plot1_single_height > 50){
        h = 50
      }else if(input$CPTAC_plot1_single_file != "pdf" & input$CPTAC_plot1_single_height <= 3000){
        h = input$CPTAC_plot1_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_plot1_single_file == "pdf" & input$CPTAC_plot1_single_width <= 50){
        w = input$CPTAC_plot1_single_width
      }else if(input$CPTAC_plot1_single_file == "pdf" & input$CPTAC_plot1_single_width > 50){
        w = 50
      }else if(input$CPTAC_plot1_single_file != "pdf" & input$CPTAC_plot1_single_width <= 3000){
        w = input$CPTAC_plot1_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_plot1_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_plot1_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_plot1_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      plotmafSummary(maf = CPTAC[[cancer_type2]][["maf"]], rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
      
      
      dev.off()
    }
  )
  output$maf2_CPTAC = renderPlot({
    cancer_type2 = cancer_type2()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = oncoplot(maf = CPTAC[[cancer_type2]][["maf"]], top = 30,legendFontSize = 2,gene_mar = 10,removeNonMutated = FALSE,legend_height=6,fontSize = 1.2)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type2())
  
  output$CPTAC_plot2_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_plot2_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","plot2", '.',"pdf")
      }else if(input$CPTAC_plot2_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","plot2", '.',"png")
      }else if(input$CPTAC_plot2_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","plot2", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","plot2", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type2 = cancer_type2()
      
      if(input$CPTAC_plot2_single_res <= 300){
        r = input$CPTAC_plot2_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_plot2_single_file == "pdf" & input$CPTAC_plot2_single_height <= 50){
        h = input$CPTAC_plot2_single_height
      }else if(input$CPTAC_plot2_single_file == "pdf" & input$CPTAC_plot2_single_height > 50){
        h = 50
      }else if(input$CPTAC_plot2_single_file != "pdf" & input$CPTAC_plot2_single_height <= 3000){
        h = input$CPTAC_plot2_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_plot2_single_file == "pdf" & input$CPTAC_plot2_single_width <= 50){
        w = input$CPTAC_plot2_single_width
      }else if(input$CPTAC_plot2_single_file == "pdf" & input$CPTAC_plot2_single_width > 50){
        w = 50
      }else if(input$CPTAC_plot2_single_file != "pdf" & input$CPTAC_plot2_single_width <= 3000){
        w = input$CPTAC_plot2_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_plot2_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_plot2_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_plot2_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      oncoplot(maf = CPTAC[[cancer_type2]][["maf"]], top = 30,legendFontSize = 2,gene_mar = 10,removeNonMutated = FALSE,legend_height=6,fontSize = 1.2)
      
      
      dev.off()
    }
  )
  
  
  output$maf3_CPTAC = renderPlot({
    cancer_type2 = cancer_type2()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = plotTiTv(titv(maf = CPTAC[[cancer_type2]][["maf"]], plot = FALSE, useSyn = TRUE),axisTextSize = c(1.3,1.3),baseFontSize = 1.3)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    
  }) %>% bindCache(cancer_type2())
  
  output$CPTAC_plot3_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_plot3_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","plot3", '.',"pdf")
      }else if(input$CPTAC_plot3_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","plot3", '.',"png")
      }else if(input$CPTAC_plot3_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","plot3", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","plot3", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type2 = cancer_type2()
      
      if(input$CPTAC_plot3_single_res <= 300){
        r = input$CPTAC_plot3_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_plot3_single_file == "pdf" & input$CPTAC_plot3_single_height <= 50){
        h = input$CPTAC_plot3_single_height
      }else if(input$CPTAC_plot3_single_file == "pdf" & input$CPTAC_plot3_single_height > 50){
        h = 50
      }else if(input$CPTAC_plot3_single_file != "pdf" & input$CPTAC_plot3_single_height <= 3000){
        h = input$CPTAC_plot3_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_plot3_single_file == "pdf" & input$CPTAC_plot3_single_width <= 50){
        w = input$CPTAC_plot3_single_width
      }else if(input$CPTAC_plot3_single_file == "pdf" & input$CPTAC_plot3_single_width > 50){
        w = 50
      }else if(input$CPTAC_plot3_single_file != "pdf" & input$CPTAC_plot3_single_width <= 3000){
        w = input$CPTAC_plot3_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_plot3_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_plot3_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_plot3_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      plotTiTv(titv(maf = CPTAC[[cancer_type2]][["maf"]], plot = FALSE, useSyn = TRUE),axisTextSize = c(1.3,1.3),baseFontSize = 1.3)
      
      
      dev.off()
    }
  )
  
  ##################################CPTAC_mutation################################################
  
  output$mut_cptac_single_uidown1 = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    tagList(
      div(selectInput(inputId = "CPTAC_lol_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_lol_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_lol_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_lol_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_lol_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_lol_single_file,{
    if(input$CPTAC_lol_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_lol_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_lol_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_lol_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_lol_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_lol_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_lol_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_cptac_single_uidown2 = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_cm_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_cm_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_cm_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_cm_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_cm_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_cm_single_file,{
    if(input$CPTAC_cm_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_cm_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_cm_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_cm_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_cm_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_cm_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_cm_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_cptac_single_uitabledown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      downloadButton('CPTAC_cm_single_tabdown',label = 'Download Table')
    )
  })
  
  observeEvent(input$sec3,{
    gene = gene_cptac()
    CPTAC_cohort = CPTAC_cohort()
    
    if(length(CPTAC_cohort$mut) <3 | length(CPTAC_cohort$wt) <3){
      closeAlert(session,"warning_cptac_single_id")
      createAlert(session, "warning_cptac_single", "warning_cptac_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_cptac_single_id")
    }
  })
  
  output$maf4_CPTAC = renderPlot({
    cancer_type2 = cancer_type2()
    gene = gene_cptac()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = lollipopPlot(maf =  CPTAC[[cancer_type2]][["maf"]], gene = gene,legendTxtSize = 2,labelOnlyUniqueDoamins = TRUE,showDomainLabel = F,
                       axisTextSize = c(1.5, 1.5),  titleSize = c(1.5, 1.2),pointSize = 2)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type2(),gene_cptac())
  
  output$CPTAC_lol_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_lol_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","lol", '.',"pdf")
      }else if(input$CPTAC_lol_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","lol", '.',"png")
      }else if(input$CPTAC_lol_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","lol", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","lol", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type2 = cancer_type2()
      gene = gene_cptac()
      
      if(input$CPTAC_lol_single_res <= 300){
        r = input$CPTAC_lol_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_lol_single_file == "pdf" & input$CPTAC_lol_single_height <= 50){
        h = input$CPTAC_lol_single_height
      }else if(input$CPTAC_lol_single_file == "pdf" & input$CPTAC_lol_single_height > 50){
        h = 50
      }else if(input$CPTAC_lol_single_file != "pdf" & input$CPTAC_lol_single_height <= 3000){
        h = input$CPTAC_lol_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_lol_single_file == "pdf" & input$CPTAC_lol_single_width <= 50){
        w = input$CPTAC_lol_single_width
      }else if(input$CPTAC_lol_single_file == "pdf" & input$CPTAC_lol_single_width > 50){
        w = 50
      }else if(input$CPTAC_lol_single_file != "pdf" & input$CPTAC_lol_single_width <= 3000){
        w = input$CPTAC_lol_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_lol_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_lol_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_lol_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      lollipopPlot(maf =  CPTAC[[cancer_type2]][["maf"]], gene = gene,legendTxtSize = 2,labelOnlyUniqueDoamins = TRUE,showDomainLabel = F,
                   axisTextSize = c(1.5, 1.5),  titleSize = c(1.5, 1.2),pointSize = 2)
      
      
      dev.off()
    }
  )
  
  
  Mut_res_CPTAC =reactive({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = Mutational_Landscape_cptac(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,mut = CPTAC_cohort$mut,wt = CPTAC_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(cancer_type2(),gene_cptac(),CPTAC_cohort())
  
  output$maf5_CPTAC = renderPlot({
    cancer_type2 = cancer_type2()
    gene = gene_cptac()
    Mut_res_CPTAC = Mut_res_CPTAC()
    
    num = length(Mut_res_CPTAC$fvsm$results$Hugo_Symbol)
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
    tmp = oncoplot(maf = Mut_res_CPTAC$tmp_maf,
                   legendFontSize = 2,
                   gene_mar = 8,
                   removeNonMutated = FALSE,
                   legend_height=6,
                   fontSize = 1.2,
                   annotationFontSize = 2,
                   annotationColor = tl,
                   annotationDat = Mut_res_CPTAC$meta,
                   clinicalFeatures = gene,
                   sortByAnnotation = T,
                   genes = Mut_res_CPTAC$fvsm$results$Hugo_Symbol[1:num],
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
    
  }) %>% bindCache(cancer_type2(),gene_cptac(),Mut_res_CPTAC())
  
  output$CPTAC_cm_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_cm_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","cm", '.',"pdf")
      }else if(input$CPTAC_cm_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","cm", '.',"png")
      }else if(input$CPTAC_cm_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","cm", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","cm", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type2 = cancer_type2()
      gene = gene_cptac()
      Mut_res_CPTAC = Mut_res_CPTAC()
      
      if(input$CPTAC_cm_single_res <= 300){
        r = input$CPTAC_cm_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_cm_single_file == "pdf" & input$CPTAC_cm_single_height <= 50){
        h = input$CPTAC_cm_single_height
      }else if(input$CPTAC_cm_single_file == "pdf" & input$CPTAC_cm_single_height > 50){
        h = 50
      }else if(input$CPTAC_cm_single_file != "pdf" & input$CPTAC_cm_single_height <= 3000){
        h = input$CPTAC_cm_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_cm_single_file == "pdf" & input$CPTAC_cm_single_width <= 50){
        w = input$CPTAC_cm_single_width
      }else if(input$CPTAC_cm_single_file == "pdf" & input$CPTAC_cm_single_width > 50){
        w = 50
      }else if(input$CPTAC_cm_single_file != "pdf" & input$CPTAC_cm_single_width <= 3000){
        w = input$CPTAC_cm_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_cm_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_cm_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_cm_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      num = length(Mut_res_CPTAC$fvsm$results$Hugo_Symbol)
      if(num > 30){
        num =30
      }
      tl = list(Mcolor)
      names(tl) = gene
      oncoplot(maf = Mut_res_CPTAC$tmp_maf,
               legendFontSize = 2,
               gene_mar = 8,
               removeNonMutated = FALSE,
               legend_height=6,
               fontSize = 1.2,
               annotationFontSize = 2,
               annotationColor = tl,
               annotationDat = Mut_res_CPTAC$meta,
               clinicalFeatures = gene,
               sortByAnnotation = T,
               genes = Mut_res_CPTAC$fvsm$results$Hugo_Symbol[1:num],
               keepGeneOrder = T)
      
      dev.off()
    }
  )
  
  
  output$compare_mutation_CPTAC = renderReactable({
    Mut_res_CPTAC = Mut_res_CPTAC()$fvsm$results
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    Mut_res_CPTAC$Hugo_Symbol = as.factor(Mut_res_CPTAC$Hugo_Symbol)
    tmp = reactable( Mut_res_CPTAC,
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
    
    
  }) %>% bindCache(Mut_res_CPTAC())
  
  output$CPTAC_cm_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      paste0(gene,"_",input$Cancer_type2,".","cm", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = Mut_res_CPTAC()$fvsm$results,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  ###################################################CPTAC immune_infiltration ################################################
  
  output$infr_cptac_single_uidown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_infr_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_infr_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infr_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infr_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_infr_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_infr_single_file,{
    if(input$CPTAC_infr_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_infr_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_infr_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_infr_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_infr_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infr_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infr_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$infp_cptac_single_uidown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_infp_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_infp_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infp_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_infp_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_infp_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_infp_single_file,{
    if(input$CPTAC_infp_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_infp_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_infp_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_infp_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_infp_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infp_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_infp_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec3,{
    gene = gene_cptac()
    CPTAC_cohort = CPTAC_cohort()
    
    if(length(CPTAC_cohort$mut) <3 | length(CPTAC_cohort$wt) <3){
      closeAlert(session,"warning2_cptac_single_id")
      createAlert(session, "warning2_cptac_single", "warning2_cptac_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning2_cptac_single_id")
    }
  })
  
  
  output$immune_infiltration1_CPTAC_rna = renderPlotly({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_infiltration_CPTAC_rna(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    tmp = ggplotly(p) %>% 
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0)) %>% 
            rangeslider(start = 0,end = 15.5)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort())
  
  output$CPTAC_infr_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_infr_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","infr", '.',"pdf")
      }else if(input$CPTAC_infr_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","infr", '.',"png")
      }else if(input$CPTAC_infr_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","infr", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","infr", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_cptac()
      cancer_type2 = cancer_type2()
      CPTAC_cohort = CPTAC_cohort()
      
      if(input$CPTAC_infr_single_res <= 300){
        r = input$CPTAC_infr_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_infr_single_file == "pdf" & input$CPTAC_infr_single_height <= 50){
        h = input$CPTAC_infr_single_height
      }else if(input$CPTAC_infr_single_file == "pdf" & input$CPTAC_infr_single_height > 50){
        h = 50
      }else if(input$CPTAC_infr_single_file != "pdf" & input$CPTAC_infr_single_height <= 3000){
        h = input$CPTAC_infr_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_infr_single_file == "pdf" & input$CPTAC_infr_single_width <= 50){
        w = input$CPTAC_infr_single_width
      }else if(input$CPTAC_infr_single_file == "pdf" & input$CPTAC_infr_single_width > 50){
        w = 50
      }else if(input$CPTAC_infr_single_file != "pdf" & input$CPTAC_infr_single_width <= 3000){
        w = input$CPTAC_infr_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_infr_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_infr_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_infr_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      print(   immune_infiltration_CPTAC_rna(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)   )
      
      
      dev.off()
    }
  )
  
  output$immune_infiltration1_CPTAC_protein = renderPlotly({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_infiltration_CPTAC_protein(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    tmp = ggplotly(p) %>% 
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0)) %>% 
            rangeslider(start = 0,end = 15.5)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort())
  
  
  output$CPTAC_infp_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_infp_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","infp", '.',"pdf")
      }else if(input$CPTAC_infp_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","infp", '.',"png")
      }else if(input$CPTAC_infp_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","infp", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","infp", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_cptac()
      cancer_type2 = cancer_type2()
      CPTAC_cohort = CPTAC_cohort()
      
      if(input$CPTAC_infp_single_res <= 300){
        r = input$CPTAC_infp_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_infp_single_file == "pdf" & input$CPTAC_infp_single_height <= 50){
        h = input$CPTAC_infp_single_height
      }else if(input$CPTAC_infp_single_file == "pdf" & input$CPTAC_infp_single_height > 50){
        h = 50
      }else if(input$CPTAC_infp_single_file != "pdf" & input$CPTAC_infp_single_height <= 3000){
        h = input$CPTAC_infp_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_infp_single_file == "pdf" & input$CPTAC_infp_single_width <= 50){
        w = input$CPTAC_infp_single_width
      }else if(input$CPTAC_infp_single_file == "pdf" & input$CPTAC_infp_single_width > 50){
        w = 50
      }else if(input$CPTAC_infp_single_file != "pdf" & input$CPTAC_infp_single_width <= 3000){
        w = input$CPTAC_infp_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_infp_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_infp_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_infp_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      print(   immune_infiltration_CPTAC_protein(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)   )
      
      
      dev.off()
    }
  )
  
  
  output$immune_infiltration2_CPTAC_rna = renderPlot({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_infiltration",selection = input$immune_cell_type_CPTAC,outlier = input$outlier2_CPTAC,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 400,height = 500) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort(),input$immune_cell_type_CPTAC,input$outlier2_CPTAC)
  
  output$immune_infiltration2_CPTAC_protein = renderPlot({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_infiltration_protein",selection = input$immune_cell_type_CPTAC,outlier = input$outlier2_CPTAC,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 400,height = 500) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort(),input$immune_cell_type_CPTAC,input$outlier2_CPTAC)
  
  ###################################################CPTAC immune_signature ################################################
  output$sigr_cptac_single_uidown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_sigr_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sigr_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigr_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigr_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sigr_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sigr_single_file,{
    if(input$CPTAC_sigr_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sigr_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sigr_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigr_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$sigp_cptac_single_uidown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_sigp_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sigp_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigp_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sigp_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sigp_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sigp_single_file,{
    if(input$CPTAC_sigp_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sigp_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sigp_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sigp_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec3,{
    gene = gene_cptac()
    CPTAC_cohort = CPTAC_cohort()
    
    if(length(CPTAC_cohort$mut) <3 | length(CPTAC_cohort$wt) <3){
      closeAlert(session,"warning3_cptac_single_id")
      createAlert(session, "warning3_cptac_single", "warning3_cptac_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning3_cptac_single_id")
    }
  })
  
  
  output$immune_signature1_CPTAC_rna = renderPlotly({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_signature_CPTAC_rna(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    tmp = ggplotly(p) %>% 
          layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0)) %>% 
          rangeslider(start = 0,end = 9)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort())
  
  output$CPTAC_sigr_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_sigr_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","sigr", '.',"pdf")
      }else if(input$CPTAC_sigr_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","sigr", '.',"png")
      }else if(input$CPTAC_sigr_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","sigr", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","sigr", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_cptac()
      cancer_type2 = cancer_type2()
      CPTAC_cohort = CPTAC_cohort()
      
      if(input$CPTAC_sigr_single_res <= 300){
        r = input$CPTAC_sigr_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_sigr_single_file == "pdf" & input$CPTAC_sigr_single_height <= 50){
        h = input$CPTAC_sigr_single_height
      }else if(input$CPTAC_sigr_single_file == "pdf" & input$CPTAC_sigr_single_height > 50){
        h = 50
      }else if(input$CPTAC_sigr_single_file != "pdf" & input$CPTAC_sigr_single_height <= 3000){
        h = input$CPTAC_sigr_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sigr_single_file == "pdf" & input$CPTAC_sigr_single_width <= 50){
        w = input$CPTAC_sigr_single_width
      }else if(input$CPTAC_sigr_single_file == "pdf" & input$CPTAC_sigr_single_width > 50){
        w = 50
      }else if(input$CPTAC_sigr_single_file != "pdf" & input$CPTAC_sigr_single_width <= 3000){
        w = input$CPTAC_sigr_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sigr_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sigr_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sigr_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      print(   immune_signature_CPTAC_rna(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)   )
      
      
      dev.off()
    }
  )
  
  
  output$immune_signature1_CPTAC_protein = renderPlotly({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_signature_CPTAC_protein(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    tmp = ggplotly(p) %>% 
            layout(boxmode = "group",legend = list(orientation = "h",x = 0.5,xanchor = "center",y = 0)) %>% 
            rangeslider(start = 0,end = 9)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort())
  
  
  output$CPTAC_sigp_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_sigp_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","sigp", '.',"pdf")
      }else if(input$CPTAC_sigp_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","sigp", '.',"png")
      }else if(input$CPTAC_sigp_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","sigp", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","sigp", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_cptac()
      cancer_type2 = cancer_type2()
      CPTAC_cohort = CPTAC_cohort()
      
      if(input$CPTAC_sigp_single_res <= 300){
        r = input$CPTAC_sigp_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_sigp_single_file == "pdf" & input$CPTAC_sigp_single_height <= 50){
        h = input$CPTAC_sigp_single_height
      }else if(input$CPTAC_sigp_single_file == "pdf" & input$CPTAC_sigp_single_height > 50){
        h = 50
      }else if(input$CPTAC_sigp_single_file != "pdf" & input$CPTAC_sigp_single_height <= 3000){
        h = input$CPTAC_sigp_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sigp_single_file == "pdf" & input$CPTAC_sigp_single_width <= 50){
        w = input$CPTAC_sigp_single_width
      }else if(input$CPTAC_sigp_single_file == "pdf" & input$CPTAC_sigp_single_width > 50){
        w = 50
      }else if(input$CPTAC_sigp_single_file != "pdf" & input$CPTAC_sigp_single_width <= 3000){
        w = input$CPTAC_sigp_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sigp_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sigp_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sigp_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(   immune_signature_CPTAC_protein(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)   )
      
      
      dev.off()
    }
  )
  
  
  output$immune_signature2_CPTAC_rna = renderPlot({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_pathway",selection = input$immune_signature_CPTAC,outlier = input$outlier3_CPTAC,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
   
   shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
   shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
   shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
   shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
   shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
   shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
   shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
   return(tmp)
   
   },width = 400,height = 500) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort(),input$immune_signature_CPTAC,input$outlier3_CPTAC)
  
  output$immune_signature2_CPTAC_protein = renderPlot({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    validate(need(input$immune_signature_CPTAC %in% rownames(CPTAC[[cancer_type2]][["immune_pathway_protein"]]),paste("The pathway",input$immune_signature_CPTAC,"does not exist in",cancer_type2,sep = " ")))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one_CPTAC(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,immune_module = "immune_pathway_protein",selection = input$immune_signature_CPTAC,outlier = input$outlier3_CPTAC,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 400,height = 500) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort(),input$immune_signature_CPTAC,input$outlier3_CPTAC)
  
  ###################################################CPTAC DEG ################################################
  
  output$DEGr_cptac_single_uitabledown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_diffr_single_tabdown',label = 'Download Table')
  })
  
  output$DEGp_cptac_single_uitabledown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_diffp_single_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec3,{
    gene = gene_cptac()
    CPTAC_cohort = CPTAC_cohort()
    
    if(length(CPTAC_cohort$mut) <3 | length(CPTAC_cohort$wt) <3){
      closeAlert(session,"warning4_cptac_single_id")
      createAlert(session, "warning4_cptac_single", "warning4_cptac_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning4_cptac_single_id")
    }
  })
  
  DEG_table_CPTAC_rna = reactive({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    min.pct_CPTAC_rna = min.pct_CPTAC_rna()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = DEG_CPTAC_rna(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,min.pct = min.pct_CPTAC_rna,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    }) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort(),min.pct_CPTAC_rna())
  
  output$DEG_tab_CPTAC_rna = renderReactable({
    DEG_table_CPTAC_rna = DEG_table_CPTAC_rna()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_CPTAC_rna,3),
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
    
  }) %>% bindCache(DEG_table_CPTAC_rna())
  
  
  output$CPTAC_diffr_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      paste0(gene,"_",input$Cancer_type2,".","diffr", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = DEG_table_CPTAC_rna(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  output$volcano_CPTAC_rna = renderPlotly({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    tab = DEG_table_CPTAC_rna()
    FC_CPTAC_rna = FC_CPTAC_rna()
    pvalue_CPTAC_rna = pvalue_CPTAC_rna()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_CPTAC_rna,pvalue = pvalue_CPTAC_rna)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(gene_cptac(),cancer_type2(),DEG_table_CPTAC_rna(),FC_CPTAC_rna(),pvalue_CPTAC_rna())
  
  #protein
  DEG_table_CPTAC_protein = reactive({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    min.pct_CPTAC_protein = min.pct_CPTAC_protein()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = DEG_CPTAC_protein(TCGA = CPTAC,gene = gene,cancer_type = cancer_type2,min.pct = min.pct_CPTAC_protein,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    }) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort(),min.pct_CPTAC_protein())
  
  output$DEG_tab_CPTAC_protein = renderReactable({
    DEG_table_CPTAC_protein = DEG_table_CPTAC_protein()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = reactable( round(DEG_table_CPTAC_protein,3),
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
  }) %>% bindCache(DEG_table_CPTAC_protein())
  
  
  output$CPTAC_diffp_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      paste0(gene,"_",input$Cancer_type2,".","diffp", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = DEG_table_CPTAC_protein(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  output$volcano_CPTAC_protein = renderPlotly({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    tab = DEG_table_CPTAC_protein()
    FC_CPTAC_protein = FC_CPTAC_protein()
    pvalue_CPTAC_protein = pvalue_CPTAC_protein()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC_CPTAC_protein,pvalue = pvalue_CPTAC_protein)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    
  }) %>% bindCache(gene_cptac(),cancer_type2(),DEG_table_CPTAC_protein(),FC_CPTAC_protein(),pvalue_CPTAC_protein())
  
  #################################################CPTAC GSEA######################################################
  
  output$GSEAr_cptac_single_uitabledown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_gsear_single_tabdown',label = 'Download Table')
  })
  
  output$GSEAp_cptac_single_uitabledown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    downloadButton('CPTAC_gseap_single_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec3,{
    gene = gene_cptac()
    CPTAC_cohort = CPTAC_cohort()
    
    if(length(CPTAC_cohort$mut) <3 | length(CPTAC_cohort$wt) <3){
      closeAlert(session,"warning5_cptac_single_id")
      createAlert(session, "warning5_cptac_single", "warning5_cptac_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning5_cptac_single_id")
    }
  })
  
  
  status_file_rna <- tempfile()
  write("", status_file_rna)
  status_file_protein <- tempfile()
  write("", status_file_protein)
  
  tmp_gsea_cptac_rna = reactive({
    diff = DEG_table_CPTAC_rna()
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
    tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_cptac_rna]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
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
  output$GSEA_tab_cptac_rna = renderReactable({
    if(tmp_gsea_cptac_rna()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader_rna() == 0){
        
        waiter_show(id = "cptac_single_rna_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader_rna(one_loader_rna()+1)
        return(NULL)
        
        
      }else if(one_loader_rna() == 1){
        
        shinyjs::disable(id = "CPTAC_gsear_single_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "cptac_single_rna_GSEA_loader")
      
      tmp_gsea_cptac_rna = tmp_gsea_cptac_rna()$get_result()@result
      tmp_gsea_cptac_rna$ID = as.factor(tmp_gsea_cptac_rna$ID)
      tmp_gsea_cptac_rna$Description = as.factor(tmp_gsea_cptac_rna$Description)
      tmp_gsea_cptac_rna$Plot = NA
      tmp_gsea_cptac_rna = tmp_gsea_cptac_rna[,c(ncol(tmp_gsea_cptac_rna),3:(ncol(tmp_gsea_cptac_rna)-2))]
      tmp = reactable( tmp_gsea_cptac_rna,defaultPageSize = 10,
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
                           Shiny.setInputValue('cptac_single_rna_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "CPTAC_gsear_single_tabdown")
      one_loader_rna(0)
      write("", status_file_rna)
      return(tmp)
    }

    
  })
  
  observeEvent(input$cptac_single_rna_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "CPTAC_gsear_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "CPTAC_gsear_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gsear_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gsear_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "CPTAC_gsear_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_cptac_rna"),type = 1)
    ))
  })
  
  observeEvent(input$CPTAC_gsear_single_file,{
    if(input$CPTAC_gsear_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_gsear_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_gsear_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gsear_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$CPTAC_gsear_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      paste0(gene,"_",input$Cancer_type2,".","gsear", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = tmp_gsea_cptac_rna()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  row_num_cptac_rna = eventReactive(input$cptac_single_rna_GSEA_show,{input$cptac_single_rna_GSEA_show$index})
  
  output$GSEA_plot_cptac_rna = renderPlot({
    my_gseaplot2(tmp_gsea_cptac_rna()$get_result(),geneSetID = tmp_gsea_cptac_rna()$get_result()@result$ID[row_num_cptac_rna()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_cptac_rna()$get_result()@result$ID[row_num_cptac_rna()])
  })
  
  
  output$CPTAC_gsear_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_gsear_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","gsear", '.',"pdf")
      }else if(input$CPTAC_gsear_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","gsear", '.',"png")
      }else if(input$CPTAC_gsear_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","gsear", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","gsear", '.',"tiff")
      }
      
    },
    content = function(file){
      
      if(input$CPTAC_gsear_single_res <= 300){
        r = input$CPTAC_gsear_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_gsear_single_file == "pdf" & input$CPTAC_gsear_single_height <= 50){
        h = input$CPTAC_gsear_single_height
      }else if(input$CPTAC_gsear_single_file == "pdf" & input$CPTAC_gsear_single_height > 50){
        h = 50
      }else if(input$CPTAC_gsear_single_file != "pdf" & input$CPTAC_gsear_single_height <= 3000){
        h = input$CPTAC_gsear_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_gsear_single_file == "pdf" & input$CPTAC_gsear_single_width <= 50){
        w = input$CPTAC_gsear_single_width
      }else if(input$CPTAC_gsear_single_file == "pdf" & input$CPTAC_gsear_single_width > 50){
        w = 50
      }else if(input$CPTAC_gsear_single_file != "pdf" & input$CPTAC_gsear_single_width <= 3000){
        w = input$CPTAC_gsear_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_gsear_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_gsear_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_gsear_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(
        my_gseaplot2(tmp_gsea_cptac_rna()$get_result(),geneSetID = tmp_gsea_cptac_rna()$get_result()@result$ID[row_num_cptac_rna()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_cptac_rna()$get_result()@result$ID[row_num_cptac_rna()])
      )
      
      
      dev.off()
    }
  )
  
  
  tmp_gsea_cptac_protein = reactive({
    diff = DEG_table_CPTAC_protein()
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
    tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_cptac_protein]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
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
  output$GSEA_tab_cptac_protein = renderReactable({
    
    if(tmp_gsea_cptac_protein()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader_protein() == 0){
        
        waiter_show(id = "cptac_single_protein_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader_protein(one_loader_protein()+1)
        return(NULL)
        
        
      }else if(one_loader_protein() == 1){
        
        shinyjs::disable(id = "CPTAC_gseap_single_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "cptac_single_protein_GSEA_loader")
      
      tmp_gsea_cptac_protein = tmp_gsea_cptac_protein()$get_result()@result
      tmp_gsea_cptac_protein$ID = as.factor(tmp_gsea_cptac_protein$ID)
      tmp_gsea_cptac_protein$Description = as.factor(tmp_gsea_cptac_protein$Description)
      tmp_gsea_cptac_protein$Plot = NA
      tmp_gsea_cptac_protein = tmp_gsea_cptac_protein[,c(ncol(tmp_gsea_cptac_protein),3:(ncol(tmp_gsea_cptac_protein)-2))]
      tmp = reactable( tmp_gsea_cptac_protein,defaultPageSize = 10,
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
                           Shiny.setInputValue('cptac_single_protein_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "CPTAC_gseap_single_tabdown")
      one_loader_protein(0)
      write("", status_file_protein)
      return(tmp)
    }

    
  })
  
  
  observeEvent(input$cptac_single_protein_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
        div(selectInput(inputId = "CPTAC_gseap_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
        div(numericInput(inputId = "CPTAC_gseap_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gseap_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(numericInput(inputId = "CPTAC_gseap_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
        div(downloadButton(outputId = "CPTAC_gseap_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_cptac_protein"),type = 1)
    ))
  })
  
  observeEvent(input$CPTAC_gseap_single_file,{
    if(input$CPTAC_gseap_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_gseap_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_gseap_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_gseap_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$CPTAC_gseap_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      paste0(gene,"_",input$Cancer_type2,".","gseap", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = tmp_gsea_cptac_protein()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  row_num_cptac_protein = eventReactive(input$cptac_single_protein_GSEA_show,{input$cptac_single_protein_GSEA_show$index})
  
  output$GSEA_plot_cptac_protein = renderPlot({
    my_gseaplot2(tmp_gsea_cptac_protein()$get_result(),geneSetID = tmp_gsea_cptac_protein()$get_result()@result$ID[row_num_cptac_protein()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_cptac_protein()$get_result()@result$ID[row_num_cptac_protein()])
  })
  
  output$CPTAC_gseap_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_gseap_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","gseap", '.',"pdf")
      }else if(input$CPTAC_gseap_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","gseap", '.',"png")
      }else if(input$CPTAC_gseap_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","gseap", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","gseap", '.',"tiff")
      }
      
    },
    content = function(file){
      
      if(input$CPTAC_gseap_single_res <= 300){
        r = input$CPTAC_gseap_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_gseap_single_file == "pdf" & input$CPTAC_gseap_single_height <= 50){
        h = input$CPTAC_gseap_single_height
      }else if(input$CPTAC_gseap_single_file == "pdf" & input$CPTAC_gseap_single_height > 50){
        h = 50
      }else if(input$CPTAC_gseap_single_file != "pdf" & input$CPTAC_gseap_single_height <= 3000){
        h = input$CPTAC_gseap_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_gseap_single_file == "pdf" & input$CPTAC_gseap_single_width <= 50){
        w = input$CPTAC_gseap_single_width
      }else if(input$CPTAC_gseap_single_file == "pdf" & input$CPTAC_gseap_single_width > 50){
        w = 50
      }else if(input$CPTAC_gseap_single_file != "pdf" & input$CPTAC_gseap_single_width <= 3000){
        w = input$CPTAC_gseap_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_gseap_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_gseap_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_gseap_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(
        
        my_gseaplot2(tmp_gsea_cptac_protein()$get_result(),geneSetID = tmp_gsea_cptac_protein()$get_result()@result$ID[row_num_cptac_protein()],
                     self.Description = "",
                     color="firebrick",
                     pvalue_table = T,
                     rel_heights = c(1, .2, 0.4),
                     title = tmp_gsea_cptac_protein()$get_result()@result$ID[row_num_cptac_protein()])
        
      )
      
      
      dev.off()
    }
  )
  
  ###################CPTAC survival###################
  
  output$sur_cptac_single_uidown = renderUI({
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "CPTAC_sur_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "CPTAC_sur_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sur_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "CPTAC_sur_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "CPTAC_sur_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$CPTAC_sur_single_file,{
    if(input$CPTAC_sur_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "CPTAC_sur_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "CPTAC_sur_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "CPTAC_sur_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "CPTAC_sur_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sur_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "CPTAC_sur_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec3,{
    gene = gene_cptac()
    CPTAC_cohort = CPTAC_cohort()
    
    if(length(CPTAC_cohort$mut) <3 | length(CPTAC_cohort$wt) <3){
      closeAlert(session,"warning6_cptac_single_id")
      createAlert(session, "warning6_cptac_single", "warning6_cptac_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning6_cptac_single_id")
    }
  })
  
  
  output$CPTAC_survival = renderPlot({
    gene = gene_cptac()
    cancer_type2 = cancer_type2()
    CPTAC_cohort = CPTAC_cohort()
    validate(need(length(CPTAC_cohort$mut) >=3 & length(CPTAC_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = CPTAC_survival(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  },width = 600,height = 600) %>% bindCache(gene_cptac(),cancer_type2(),CPTAC_cohort())
  
  output$CPTAC_sur_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_cptac()
      
      if(input$CPTAC_sur_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type2,".","sur", '.',"pdf")
      }else if(input$CPTAC_sur_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type2,".","sur", '.',"png")
      }else if(input$CPTAC_sur_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type2,".","sur", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type2,".","sur", '.',"tiff")
      }
      
    },
    content = function(file){
      
      gene = gene_cptac()
      cancer_type2 = cancer_type2()
      CPTAC_cohort = CPTAC_cohort()
      
      if(input$CPTAC_sur_single_res <= 300){
        r = input$CPTAC_sur_single_res
      }else{
        r = 300
      }
      if(input$CPTAC_sur_single_file == "pdf" & input$CPTAC_sur_single_height <= 50){
        h = input$CPTAC_sur_single_height
      }else if(input$CPTAC_sur_single_file == "pdf" & input$CPTAC_sur_single_height > 50){
        h = 50
      }else if(input$CPTAC_sur_single_file != "pdf" & input$CPTAC_sur_single_height <= 3000){
        h = input$CPTAC_sur_single_height
      }else{
        h = 3000
      }
      if(input$CPTAC_sur_single_file == "pdf" & input$CPTAC_sur_single_width <= 50){
        w = input$CPTAC_sur_single_width
      }else if(input$CPTAC_sur_single_file == "pdf" & input$CPTAC_sur_single_width > 50){
        w = 50
      }else if(input$CPTAC_sur_single_file != "pdf" & input$CPTAC_sur_single_width <= 3000){
        w = input$CPTAC_sur_single_width
      }else{
        w = 3000
      }
      
      if(input$CPTAC_sur_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$CPTAC_sur_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$CPTAC_sur_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      print(    CPTAC_survival(TCGA = CPTAC,cancer_type = cancer_type2,gene = gene,mut_patient_id = CPTAC_cohort$mut,wt_patient_id = CPTAC_cohort$wt)   )
      
      
      dev.off()
    }
  )
  
  
}