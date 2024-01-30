TCGA_server <- function(input,output,session,gene_tcga,cancer_type,TCGA,pathway_database,Mut_type_tcga_single,Wild_type_tcga_single,min.pct,FC,pvalue){
  
  
  TCGA_cohort = reactive({
    cancer_type = cancer_type()
    gene = gene_tcga()
    Mut_type_tcga_single = Mut_type_tcga_single()
    Wild_type_tcga_single = Wild_type_tcga_single()
    
    TCGA_cohort_cal(TCGA = TCGA,cancer_type = cancer_type,gene = gene,Mut_type = Mut_type_tcga_single,Wild_type = Wild_type_tcga_single)
  }) %>% bindCache(cancer_type(),gene_tcga(),Mut_type_tcga_single(),Wild_type_tcga_single())
  
  ###################################TCGA_overview################################################
  
  output$over_tcga_single_uidown1 = renderUI({
    cancer_type = cancer_type()
    tagList(
      div(selectInput(inputId = "TCGA_plot1_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_plot1_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_plot1_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_plot1_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_plot1_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      )
  })
  
  observeEvent(input$TCGA_plot1_single_file,{
    if(input$TCGA_plot1_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_plot1_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_plot1_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_plot1_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_plot1_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_plot1_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_plot1_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$over_tcga_single_uidown2 = renderUI({
  cancer_type = cancer_type()
  tagList(
    div(selectInput(inputId = "TCGA_plot2_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
    div(numericInput(inputId = "TCGA_plot2_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
    div(numericInput(inputId = "TCGA_plot2_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
    div(numericInput(inputId = "TCGA_plot2_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
    div(downloadButton(outputId = "TCGA_plot2_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
})
  
  observeEvent(input$TCGA_plot2_single_file,{
    if(input$TCGA_plot2_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_plot2_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_plot2_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_plot2_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_plot2_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_plot2_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_plot2_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })

  output$over_tcga_single_uidown3 = renderUI({
  cancer_type = cancer_type()
  tagList(
    div(selectInput(inputId = "TCGA_plot3_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
    div(numericInput(inputId = "TCGA_plot3_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
    div(numericInput(inputId = "TCGA_plot3_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
    div(numericInput(inputId = "TCGA_plot3_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
    div(downloadButton(outputId = "TCGA_plot3_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
  )
})
  
  observeEvent(input$TCGA_plot3_single_file,{
    if(input$TCGA_plot3_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_plot3_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_plot3_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_plot3_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_plot3_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_plot3_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_plot3_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  
  output$maf1 = renderPlot({
    cancer_type = cancer_type()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = plotmafSummary(maf = TCGA[[cancer_type]][["maf"]], rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type())
  
  output$TCGA_plot1_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_plot1_single_file == 'pdf'){
            paste0(gene,"_",input$Cancer_type,".","plot1", '.',"pdf")
      }else if(input$TCGA_plot1_single_file == 'png'){
            paste0(gene,"_",input$Cancer_type,".","plot1", '.',"png")
      }else if(input$TCGA_plot1_single_file == 'jpeg'){
            paste0(gene,"_",input$Cancer_type,".","plot1", '.',"jpeg")
      }else{
            paste0(gene,"_",input$Cancer_type,".","plot1", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type = cancer_type()
      
      if(input$TCGA_plot1_single_res <= 300){
        r = input$TCGA_plot1_single_res
      }else{
        r = 300
      }
      if(input$TCGA_plot1_single_file == "pdf" & input$TCGA_plot1_single_height <= 50){
        h = input$TCGA_plot1_single_height
      }else if(input$TCGA_plot1_single_file == "pdf" & input$TCGA_plot1_single_height > 50){
        h = 50
      }else if(input$TCGA_plot1_single_file != "pdf" & input$TCGA_plot1_single_height <= 3000){
        h = input$TCGA_plot1_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_plot1_single_file == "pdf" & input$TCGA_plot1_single_width <= 50){
        w = input$TCGA_plot1_single_width
      }else if(input$TCGA_plot1_single_file == "pdf" & input$TCGA_plot1_single_width > 50){
        w = 50
      }else if(input$TCGA_plot1_single_file != "pdf" & input$TCGA_plot1_single_width <= 3000){
        w = input$TCGA_plot1_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_plot1_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_plot1_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_plot1_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      plotmafSummary(maf = TCGA[[cancer_type]][["maf"]], rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
      
      
      dev.off()
    }
  )
  
  output$maf2 = renderPlot({
    cancer_type = cancer_type()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = oncoplot(maf = TCGA[[cancer_type]][["maf"]],top = 30,legendFontSize = 2,gene_mar = 10,removeNonMutated = FALSE,legend_height=6,fontSize = 1.2)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type())
  
  output$TCGA_plot2_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_plot2_single_file == 'pdf'){
            paste0(gene,"_",input$Cancer_type,".","plot2", '.',"pdf")
      }else if(input$TCGA_plot2_single_file == 'png'){
            paste0(gene,"_",input$Cancer_type,".","plot2", '.',"png")
      }else if(input$TCGA_plot2_single_file == 'jpeg'){
            paste0(gene,"_",input$Cancer_type,".","plot2", '.',"jpeg")
      }else{
            paste0(gene,"_",input$Cancer_type,".","plot2", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type = cancer_type()
      
      if(input$TCGA_plot2_single_res <= 300){
        r = input$TCGA_plot2_single_res
      }else{
        r = 300
      }
      if(input$TCGA_plot2_single_file == "pdf" & input$TCGA_plot2_single_height <= 50){
        h = input$TCGA_plot2_single_height
      }else if(input$TCGA_plot2_single_file == "pdf" & input$TCGA_plot2_single_height > 50){
        h = 50
      }else if(input$TCGA_plot2_single_file != "pdf" & input$TCGA_plot2_single_height <= 3000){
        h = input$TCGA_plot2_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_plot2_single_file == "pdf" & input$TCGA_plot2_single_width <= 50){
        w = input$TCGA_plot2_single_width
      }else if(input$TCGA_plot2_single_file == "pdf" & input$TCGA_plot2_single_width > 50){
        w = 50
      }else if(input$TCGA_plot2_single_file != "pdf" & input$TCGA_plot2_single_width <= 3000){
        w = input$TCGA_plot2_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_plot2_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_plot2_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_plot2_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }

      
      oncoplot(maf = TCGA[[cancer_type]][["maf"]],top = 30,legendFontSize =2,gene_mar = 10,removeNonMutated = FALSE,legend_height=6,fontSize = 1.2)
      
      
      dev.off()
    }
  )
  
  output$maf3 = renderPlot({
    cancer_type = cancer_type()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = plotTiTv(titv(maf = TCGA[[cancer_type]][["maf"]], plot = FALSE, useSyn = TRUE),axisTextSize = c(1.3,1.3),baseFontSize = 1.3)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type())
  
  output$TCGA_plot3_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_plot3_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","plot3", '.',"pdf")
      }else if(input$TCGA_plot3_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type,".","plot3", '.',"png")
      }else if(input$TCGA_plot3_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","plot3", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","plot3", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type = cancer_type()
      
      if(input$TCGA_plot3_single_res <= 300){
        r = input$TCGA_plot3_single_res
      }else{
        r = 300
      }
      if(input$TCGA_plot3_single_file == "pdf" & input$TCGA_plot3_single_height <= 50){
        h = input$TCGA_plot3_single_height
      }else if(input$TCGA_plot3_single_file == "pdf" & input$TCGA_plot3_single_height > 50){
        h = 50
      }else if(input$TCGA_plot3_single_file != "pdf" & input$TCGA_plot3_single_height <= 3000){
        h = input$TCGA_plot3_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_plot3_single_file == "pdf" & input$TCGA_plot3_single_width <= 50){
        w = input$TCGA_plot3_single_width
      }else if(input$TCGA_plot3_single_file == "pdf" & input$TCGA_plot3_single_width > 50){
        w = 50
      }else if(input$TCGA_plot3_single_file != "pdf" & input$TCGA_plot3_single_width <= 3000){
        w = input$TCGA_plot3_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_plot3_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_plot3_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_plot3_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      plotTiTv(titv(maf = TCGA[[cancer_type]][["maf"]], plot = FALSE, useSyn = TRUE),axisTextSize = c(1.3,1.3),baseFontSize = 1.3)
      
      
      dev.off()
    }
  )
  
  ##################################TCGA_mutation################################################
  output$mut_tcga_single_uidown1 = renderUI({
    TCGA_cohort = TCGA_cohort()
    tagList(
      div(selectInput(inputId = "TCGA_lol_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_lol_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_lol_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_lol_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_lol_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_lol_single_file,{
    if(input$TCGA_lol_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_lol_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_lol_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_lol_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_lol_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_lol_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_lol_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_tcga_single_uidown2 = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_cm_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_cm_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_cm_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_cm_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_cm_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_cm_single_file,{
    if(input$TCGA_cm_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_cm_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_cm_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_cm_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_cm_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_cm_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_cm_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$mut_tcga_single_uitabledown = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    tagList(
      downloadButton('TCGA_cm_single_tabdown',label = 'Download Table')
    )
  })
  
  
  observeEvent(input$sec2,{
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    
    if(length(TCGA_cohort$mut) <3 | length(TCGA_cohort$wt) <3){
      closeAlert(session,"warning_tcga_single_id")
      createAlert(session, "warning_tcga_single", "warning_tcga_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning_tcga_single_id")
    }
  })
  
  output$maf4 = renderPlot({
    cancer_type = cancer_type()
    gene = gene_tcga()
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = lollipopPlot(maf =  TCGA[[cancer_type]][["maf"]], gene = gene,legendTxtSize = 2,labelOnlyUniqueDoamins = TRUE,showDomainLabel = F,
                       axisTextSize = c(1.5, 1.5),  titleSize = c(1.5, 1.2),pointSize = 2)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type(),gene_tcga())
  
  output$TCGA_lol_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_lol_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","lol", '.',"pdf")
      }else if(input$TCGA_lol_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type,".","lol", '.',"png")
      }else if(input$TCGA_lol_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","lol", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","lol", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type = cancer_type()
      gene = gene_tcga()
      
      if(input$TCGA_lol_single_res <= 300){
        r = input$TCGA_lol_single_res
      }else{
        r = 300
      }
      if(input$TCGA_lol_single_file == "pdf" & input$TCGA_lol_single_height <= 50){
        h = input$TCGA_lol_single_height
      }else if(input$TCGA_lol_single_file == "pdf" & input$TCGA_lol_single_height > 50){
        h = 50
      }else if(input$TCGA_lol_single_file != "pdf" & input$TCGA_lol_single_height <= 3000){
        h = input$TCGA_lol_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_lol_single_file == "pdf" & input$TCGA_lol_single_width <= 50){
        w = input$TCGA_lol_single_width
      }else if(input$TCGA_lol_single_file == "pdf" & input$TCGA_lol_single_width > 50){
        w = 50
      }else if(input$TCGA_lol_single_file != "pdf" & input$TCGA_lol_single_width <= 3000){
        w = input$TCGA_lol_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_lol_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_lol_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_lol_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      lollipopPlot(maf =  TCGA[[cancer_type]][["maf"]], gene = gene,legendTxtSize = 2,labelOnlyUniqueDoamins = TRUE,showDomainLabel = F,
                   axisTextSize = c(1.5, 1.5),  titleSize = c(1.5, 1.2),pointSize = 2)
      
      
      dev.off()
    }
  )
  
  Mut_res =reactive({
    cancer_type = cancer_type()
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = Mutational_Landscape(TCGA = TCGA,cancer_type = cancer_type,gene = gene,mut = TCGA_cohort$mut,wt = TCGA_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
  }) %>% bindCache(cancer_type(),gene_tcga(),TCGA_cohort())
  
  
  output$maf5 = renderPlot({
    cancer_type = cancer_type()
    gene = gene_tcga()
    Mut_res = Mut_res()
    tmp_maf = subsetMaf(maf = TCGA[[cancer_type]][["maf"]],tsb = Mut_res$meta$Tumor_Sample_Barcode)
    num = length(Mut_res$fvsm$results$Hugo_Symbol)
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
    tmp = oncoplot(maf = tmp_maf,
                   legendFontSize = 2,
                   gene_mar = 8,
                   removeNonMutated = FALSE,
                   legend_height=6,
                   fontSize = 1.2,
                   clinicalFeatures = gene,
                   annotationColor = tl,
                   annotationDat = Mut_res$meta,
                   annotationFontSize = 2,
                   sortByAnnotation = T,
                   genes = Mut_res$fvsm$results$Hugo_Symbol[1:num],
                   keepGeneOrder = T)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(cancer_type(),gene_tcga(), Mut_res())
  
  
  output$TCGA_cm_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_cm_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","cm", '.',"pdf")
      }else if(input$TCGA_cm_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type,".","cm", '.',"png")
      }else if(input$TCGA_cm_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","cm", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","cm", '.',"tiff")
      }
      
    },
    content = function(file){
      cancer_type = cancer_type()
      gene = gene_tcga()
      Mut_res = Mut_res()
      
      if(input$TCGA_cm_single_res <= 300){
        r = input$TCGA_cm_single_res
      }else{
        r = 300
      }
      if(input$TCGA_cm_single_file == "pdf" & input$TCGA_cm_single_height <= 50){
        h = input$TCGA_cm_single_height
      }else if(input$TCGA_cm_single_file == "pdf" & input$TCGA_cm_single_height > 50){
        h = 50
      }else if(input$TCGA_cm_single_file != "pdf" & input$TCGA_cm_single_height <= 3000){
        h = input$TCGA_cm_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_cm_single_file == "pdf" & input$TCGA_cm_single_width <= 50){
        w = input$TCGA_cm_single_width
      }else if(input$TCGA_cm_single_file == "pdf" & input$TCGA_cm_single_width > 50){
        w = 50
      }else if(input$TCGA_cm_single_file != "pdf" & input$TCGA_cm_single_width <= 3000){
        w = input$TCGA_cm_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_cm_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_cm_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_cm_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      

      tmp_maf = subsetMaf(maf = TCGA[[cancer_type]][["maf"]],tsb = Mut_res$meta$Tumor_Sample_Barcode)
      num = length(Mut_res$fvsm$results$Hugo_Symbol)
      if(num > 30){
        num =30
      }
      
      tl = list(Mcolor)
      names(tl) = gene
      tmp = oncoplot(maf = tmp_maf,
                     legendFontSize = 2,
                     gene_mar = 8,
                     removeNonMutated = FALSE,
                     legend_height=6,
                     fontSize = 1.2,
                     clinicalFeatures = gene,
                     annotationColor = tl,
                     annotationDat = Mut_res$meta,
                     annotationFontSize = 2,
                     sortByAnnotation = T,
                     genes = Mut_res$fvsm$results$Hugo_Symbol[1:num],
                     keepGeneOrder = T)      
      
      dev.off()
    }
  )
  
  output$compare_mutation = renderReactable({
    Mut_res = Mut_res()$fvsm$results
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    Mut_res$Hugo_Symbol = as.factor(Mut_res$Hugo_Symbol)
    tmp = reactable( Mut_res,
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
    
  }) %>% bindCache( Mut_res() )
  
  output$TCGA_cm_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      paste0(gene,"_",input$Cancer_type,".","cm", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = Mut_res()$fvsm$results,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  ###################################################TCGA immune_infiltration ################################################
  
  output$inf_tcga_single_uidown = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_inf_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_inf_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_inf_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_inf_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_inf_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_inf_single_file,{
    if(input$TCGA_inf_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_inf_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_inf_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_inf_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_inf_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_inf_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_inf_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec2,{
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    
    if(length(TCGA_cohort$mut) <3 | length(TCGA_cohort$wt) <3){
      closeAlert(session,"warning2_tcga_single_id")
      createAlert(session, "warning2_tcga_single", "warning2_tcga_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning2_tcga_single_id")
    }
  })
  
  output$immune_infiltration1 = renderPlotly({
    gene = gene_tcga()
    cancer_type = cancer_type()
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_infiltration(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt)
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
    
    }) %>% bindCache(gene_tcga(),cancer_type(),TCGA_cohort())
  
  output$TCGA_inf_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_inf_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","inf", '.',"pdf")
      }else if(input$TCGA_inf_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type,".","inf", '.',"png")
      }else if(input$TCGA_inf_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","inf", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","inf", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_tcga()
      cancer_type = cancer_type()
      TCGA_cohort = TCGA_cohort()
      
      if(input$TCGA_inf_single_res <= 300){
        r = input$TCGA_inf_single_res
      }else{
        r = 300
      }
      if(input$TCGA_inf_single_file == "pdf" & input$TCGA_inf_single_height <= 50){
        h = input$TCGA_inf_single_height
      }else if(input$TCGA_inf_single_file == "pdf" & input$TCGA_inf_single_height > 50){
        h = 50
      }else if(input$TCGA_inf_single_file != "pdf" & input$TCGA_inf_single_height <= 3000){
        h = input$TCGA_inf_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_inf_single_file == "pdf" & input$TCGA_inf_single_width <= 50){
        w = input$TCGA_inf_single_width
      }else if(input$TCGA_inf_single_file == "pdf" & input$TCGA_inf_single_width > 50){
        w = 50
      }else if(input$TCGA_inf_single_file != "pdf" & input$TCGA_inf_single_width <= 3000){
        w = input$TCGA_inf_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_inf_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_inf_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_inf_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      print(immune_infiltration(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt))
      
      dev.off()
    }
  )
  
  output$immune_infiltration2 = renderPlot({
    gene = gene_tcga()
    cancer_type = cancer_type()
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one(TCGA = TCGA,gene = gene,cancer_type = cancer_type,immune_module = "immune_infiltration",selection = input$immune_cell_type,outlier = input$outlier2,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 600,height = 600) %>% bindCache(gene_tcga(),cancer_type(),TCGA_cohort(),input$immune_cell_type,input$outlier2)
  ###################################################TCGA immune_signature ################################################
  
  output$sig_tcga_single_uidown = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_sig_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_sig_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sig_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sig_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_sig_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_sig_single_file,{
    if(input$TCGA_sig_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_sig_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_sig_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_sig_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_sig_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sig_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sig_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec2,{
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    
    if(length(TCGA_cohort$mut) <3 | length(TCGA_cohort$wt) <3){
      closeAlert(session,"warning3_tcga_single_id")
      createAlert(session, "warning3_tcga_single", "warning3_tcga_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning3_tcga_single_id")
    }
  })
  
  
  output$immune_signature1 = renderPlotly({
    gene = gene_tcga()
    cancer_type = cancer_type()
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    p = immune_signature(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt)
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
    
    }) %>% bindCache(gene_tcga(),cancer_type(),TCGA_cohort())
  
  
  
  
  output$TCGA_sig_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_sig_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","sig", '.',"pdf")
      }else if(input$TCGA_sig_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type,".","sig", '.',"png")
      }else if(input$TCGA_sig_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","sig", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","sig", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_tcga()
      cancer_type = cancer_type()
      TCGA_cohort = TCGA_cohort()
      
      if(input$TCGA_sig_single_res <= 300){
        r = input$TCGA_sig_single_res
      }else{
        r = 300
      }
      if(input$TCGA_sig_single_file == "pdf" & input$TCGA_sig_single_height <= 50){
        h = input$TCGA_sig_single_height
      }else if(input$TCGA_sig_single_file == "pdf" & input$TCGA_sig_single_height > 50){
        h = 50
      }else if(input$TCGA_sig_single_file != "pdf" & input$TCGA_sig_single_height <= 3000){
        h = input$TCGA_sig_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_sig_single_file == "pdf" & input$TCGA_sig_single_width <= 50){
        w = input$TCGA_sig_single_width
      }else if(input$TCGA_sig_single_file == "pdf" & input$TCGA_sig_single_width > 50){
        w = 50
      }else if(input$TCGA_sig_single_file != "pdf" & input$TCGA_sig_single_width <= 3000){
        w = input$TCGA_sig_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_sig_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_sig_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_sig_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      print(immune_signature(TCGA = TCGA,gene = gene,cancer_type = cancer_type,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt))
      
      dev.off()
    }
  )
  
  output$immune_signature2 = renderPlot({
    gene = gene_tcga()
    cancer_type = cancer_type()
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = immune_one(TCGA = TCGA,gene = gene,cancer_type = cancer_type,immune_module = "immune_pathway",selection = input$immune_signature,outlier = input$outlier3,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt)
  
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    },width = 600,height = 600) %>% bindCache(gene_tcga(),cancer_type(),TCGA_cohort(),input$immune_signature,input$outlier3)
  
  
  ###################################################TCGA DEG ################################################
  
  output$DEG_tcga_single_uitabledown = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    downloadButton('TCGA_diff_single_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec2,{
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    
    if(length(TCGA_cohort$mut) <3 | length(TCGA_cohort$wt) <3){
      closeAlert(session,"warning4_tcga_single_id")
      createAlert(session, "warning4_tcga_single", "warning4_tcga_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning4_tcga_single_id")
    }
  })
  
  DEG_table = reactive({
    gene = gene_tcga()
    cancer_type = cancer_type()
    TCGA_cohort = TCGA_cohort()
    min.pct = min.pct()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = DEG(TCGA = TCGA,gene = gene,cancer_type = cancer_type,min.pct = min.pct,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
    }) %>% bindCache(gene_tcga(),cancer_type(),TCGA_cohort(),min.pct())
  
  output$DEG_tab = renderReactable({
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
  
  output$TCGA_diff_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      paste0(gene,"_",input$Cancer_type,".","diff", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = DEG_table(),file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  output$volcano = renderPlotly({
    tab = DEG_table()
    FC = FC()
    pvalue = pvalue()
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = volcano_plot(tab = tab,FC = FC,pvalue = pvalue)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)
    
  }) %>% bindCache(DEG_table(),FC(),pvalue())
  
  #################################################TCGA GSEA######################################################
  
  output$GSEA_tcga_single_uitabledown = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    downloadButton('TCGA_gsea_single_tabdown',label = 'Download Table')
  })
  
  observeEvent(input$sec2,{
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    
    if(length(TCGA_cohort$mut) <3 | length(TCGA_cohort$wt) <3){
      closeAlert(session,"warning5_tcga_single_id")
      createAlert(session, "warning5_tcga_single", "warning5_tcga_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning5_tcga_single_id")
    }
  })
  
  
  status_file <- tempfile()
  write("", status_file)
  
  tmp_gsea_tcga = reactive({
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
    
    installr::kill_pid(pid = scan(status_file))
    tmp = r_bg(func = myGSEA,args = list(geneList = FC,TERM2GENE=pathway_database[[input$pathway_gsea_tcga]][,c(1,3)],pvalueCutoff = 0.1),supervise = TRUE)
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
  output$GSEA_tab_tcga = renderReactable({
    
    if(tmp_gsea_tcga()$is_alive()){
      
      invalidateLater(millis = 1000, session = session)
      
      if(one_loader() == 0){
        
        waiter_show(id = "tcga_single_GSEA_loader",html = tagList(spin_flower(),h4("GSEA running..."),h5("Please wait for a minute")), color = "black")
        one_loader(one_loader()+1)
        return(NULL)
        
        
      }else if(one_loader() == 1){
        
        shinyjs::disable(id = "TCGA_gsea_single_tabdown")
        
      }
      
    }else{
      
      waiter_hide(id = "tcga_single_GSEA_loader")
      
      tmp_gsea_tcga = tmp_gsea_tcga()$get_result()@result
      tmp_gsea_tcga$ID = as.factor(tmp_gsea_tcga$ID)
      tmp_gsea_tcga$Description = as.factor(tmp_gsea_tcga$Description)
      tmp_gsea_tcga$Plot = NA
      tmp_gsea_tcga = tmp_gsea_tcga[,c(ncol(tmp_gsea_tcga),3:(ncol(tmp_gsea_tcga)-2))]
      tmp = reactable( tmp_gsea_tcga,defaultPageSize = 10,
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
                           Shiny.setInputValue('tcga_single_GSEA_show', { index: rowInfo.index + 1 }, { priority: 'event' })
                         }
                       }"
                       )
                       
      )
      
      
      shinyjs::enable(id = "TCGA_gsea_single_tabdown")
      one_loader(0)
      write("", status_file)
      return(tmp)
    }
  }) 
  
  observeEvent(input$tcga_single_GSEA_show, {
    showModal(modalDialog(
      title = "GSEA Plot",
      easyClose = TRUE,
      size = "l",
      footer = tagList(
               div(selectInput(inputId = "TCGA_gsea_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
               div(numericInput(inputId = "TCGA_gsea_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
               div(numericInput(inputId = "TCGA_gsea_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
               div(numericInput(inputId = "TCGA_gsea_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
               div(downloadButton(outputId = "TCGA_gsea_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
      ),
      withSpinner(plotOutput(outputId = "GSEA_plot_tcga"),type = 1)
    ))
  })
  
  observeEvent(input$TCGA_gsea_single_file,{
    if(input$TCGA_gsea_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_gsea_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_gsea_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_gsea_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_gsea_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_gsea_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_gsea_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  output$TCGA_gsea_single_tabdown = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      paste0(gene,"_",input$Cancer_type,".","gsea", '.',"csv")
      
    },
    content = function(file){
      
      write.csv(x = tmp_gsea_tcga()$get_result()@result,file = file, sep = ',', col.names = T, row.names = T, quote = F)
      
    }
  )
  
  
  row_num_tcga = eventReactive(input$tcga_single_GSEA_show,{input$tcga_single_GSEA_show$index})
  
  output$GSEA_plot_tcga = renderPlot({
    my_gseaplot2(tmp_gsea_tcga()$get_result(),geneSetID = tmp_gsea_tcga()$get_result()@result$ID[row_num_tcga()],
                 self.Description = "",
                 color="firebrick",
                 pvalue_table = T,
                 rel_heights = c(1, .2, 0.4),
                 title = tmp_gsea_tcga()$get_result()@result$ID[row_num_tcga()])

  })
  
  output$TCGA_gsea_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      message(input$TCGA_gsea_single_file)
      if(input$TCGA_gsea_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","gsea", '.',"pdf")
      }else if(input$TCGA_gsea_single_file == 'png'){
        message(input$TCGA_gsea_single_file)
        paste0(gene,"_",input$Cancer_type,".","gsea", '.',"png")
      }else if(input$TCGA_gsea_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","gsea", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","gsea", '.',"tiff")
      }
      
    },
    content = function(file){
      
      if(input$TCGA_gsea_single_res <= 300){
        r = input$TCGA_gsea_single_res
      }else{
        r = 300
      }
      if(input$TCGA_gsea_single_file == "pdf" & input$TCGA_gsea_single_height <= 50){
        h = input$TCGA_gsea_single_height
      }else if(input$TCGA_gsea_single_file == "pdf" & input$TCGA_gsea_single_height > 50){
        h = 50
      }else if(input$TCGA_gsea_single_file != "pdf" & input$TCGA_gsea_single_height <= 3000){
        h = input$TCGA_gsea_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_gsea_single_file == "pdf" & input$TCGA_gsea_single_width <= 50){
        w = input$TCGA_gsea_single_width
      }else if(input$TCGA_gsea_single_file == "pdf" & input$TCGA_gsea_single_width > 50){
        w = 50
      }else if(input$TCGA_gsea_single_file != "pdf" & input$TCGA_gsea_single_width <= 3000){
        w = input$TCGA_gsea_single_width
      }else{
        w = 3000
      }

      
      if(input$TCGA_gsea_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_gsea_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_gsea_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      print(    my_gseaplot2(tmp_gsea_tcga()$get_result(),geneSetID = tmp_gsea_tcga()$get_result()@result$ID[row_num_tcga()],
                             self.Description = "",
                             color="firebrick",
                             pvalue_table = T,
                             rel_heights = c(1, .2, 0.4),
                             title = tmp_gsea_tcga()$get_result()@result$ID[row_num_tcga()])
                )
      
      dev.off()
    }
  )
  
  
  ###################TCGA survival###################
  
  output$sur_tcga_single_uidown = renderUI({
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    tagList(
      div(selectInput(inputId = "TCGA_sur_single_file",label = "Choose a format",choices = c("png","tiff","jpeg","pdf"),selected = "png",multiple = F,width = "100%"),style = "display:inline-block;width:20%;vertical-align:top;margin-left:20px;"),
      div(numericInput(inputId = "TCGA_sur_single_res",label = "Resolution",min = 100,max = 300,value = 100,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sur_single_width",label = "Width",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(numericInput(inputId = "TCGA_sur_single_height",label = "Height",min = 500,max = 3000,value = 1000,step = 10,width = "100%"),style = "display:inline-block;width:20%;"),
      div(downloadButton(outputId = "TCGA_sur_single_down",label = "Download Plot"),style = "display:inline-block;width:10%;vertical-align: middle;")
    )
  })
  
  observeEvent(input$TCGA_sur_single_file,{
    if(input$TCGA_sur_single_file == "pdf"){
      updateNumericInput(session = session,inputId = "TCGA_sur_single_res",value = 0,min = 0,max = 0)
      updateNumericInput(session = session,inputId = "TCGA_sur_single_width",value = 10,min = 1,max = 30,step = 1)
      updateNumericInput(session = session,inputId = "TCGA_sur_single_height",value = 10,min = 1,max = 30,step = 1)
    }else{
      updateNumericInput(session = session,inputId = "TCGA_sur_single_res",min = 100,max = 300,value = 100,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sur_single_width",min = 500,max = 3000,value = 1000,step = 10)
      updateNumericInput(session = session,inputId = "TCGA_sur_single_height",min = 500,max = 3000,value = 1000,step = 10)
    }
  })
  
  observeEvent(input$sec2,{
    gene = gene_tcga()
    TCGA_cohort = TCGA_cohort()
    
    if(length(TCGA_cohort$mut) <3 | length(TCGA_cohort$wt) <3){
      closeAlert(session,"warning6_tcga_single_id")
      createAlert(session, "warning6_tcga_single", "warning6_tcga_single_id", title = "Warning",style = "danger",
                  content = paste("The number of patients with",gene,"mutation or wildtype","<3",sep = " "), append = FALSE)
    }else{
      closeAlert(session,"warning6_tcga_single_id")
    }
  })
  
  
  width = reactive({    if(cancer_type() == "LAML"){500}else{1600}})
  output$TCGA_survival = renderPlot({
    gene = gene_tcga()
    cancer_type = cancer_type()
    TCGA_cohort = TCGA_cohort()
    validate(need(length(TCGA_cohort$mut) >=3 & length(TCGA_cohort$wt) >= 3,message = FALSE))
    shinyjs::disable("sec");shinyjs::disable("sec2");shinyjs::disable("sec3");shinyjs::disable("sec_pm");shinyjs::disable("sec2_pm");shinyjs::disable("sec3_pm");shinyjs::disable("sec2_subtype");shinyjs::disable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    
    tmp = TCGA_survival(TCGA = TCGA,cancer_type = cancer_type,gene = gene,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt)
    
    shinyjs::enable("sec");shinyjs::enable("sec2");shinyjs::enable("sec3");shinyjs::enable("sec_pm");shinyjs::enable("sec2_pm");shinyjs::enable("sec3_pm");shinyjs::enable("sec2_subtype");shinyjs::enable("sec3_subtype");
    shinyjs::removeClass(id = "button_ref_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_single",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_single",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_ref_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_tcga_pm",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_pm",class = "fa fa-spinner fa-spin");
    shinyjs::removeClass(id = "button_tcga_subtype",class = "fa fa-spinner fa-spin");shinyjs::removeClass(id = "button_cptac_subtype",class = "fa fa-spinner fa-spin");
    shinyjs::addClass(id = "button_ref_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_single",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_single",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_ref_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_tcga_pm",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_pm",class = "fa fa-arrow-alt-circle-up");
    shinyjs::addClass(id = "button_tcga_subtype",class = "fa fa-arrow-alt-circle-up");shinyjs::addClass(id = "button_cptac_subtype",class = "fa fa-arrow-alt-circle-up");
    return(tmp)

  },width = width,height = 600) %>% bindCache(gene_tcga(),cancer_type(),TCGA_cohort())
  
  output$TCGA_sur_single_down = downloadHandler(
    
    filename = function(){
      gene = gene_tcga()
      
      if(input$TCGA_sur_single_file == 'pdf'){
        paste0(gene,"_",input$Cancer_type,".","sur", '.',"pdf")
      }else if(input$TCGA_sur_single_file == 'png'){
        paste0(gene,"_",input$Cancer_type,".","sur", '.',"png")
      }else if(input$TCGA_sur_single_file == 'jpeg'){
        paste0(gene,"_",input$Cancer_type,".","sur", '.',"jpeg")
      }else{
        paste0(gene,"_",input$Cancer_type,".","sur", '.',"tiff")
      }
      
    },
    content = function(file){
      gene = gene_tcga()
      cancer_type = cancer_type()
      TCGA_cohort = TCGA_cohort()
      
      if(input$TCGA_sur_single_res <= 300){
        r = input$TCGA_sur_single_res
      }else{
        r = 300
      }
      if(input$TCGA_sur_single_file == "pdf" & input$TCGA_sur_single_height <= 50){
        h = input$TCGA_sur_single_height
      }else if(input$TCGA_sur_single_file == "pdf" & input$TCGA_sur_single_height > 50){
        h = 50
      }else if(input$TCGA_sur_single_file != "pdf" & input$TCGA_sur_single_height <= 3000){
        h = input$TCGA_sur_single_height
      }else{
        h = 3000
      }
      if(input$TCGA_sur_single_file == "pdf" & input$TCGA_sur_single_width <= 50){
        w = input$TCGA_sur_single_width
      }else if(input$TCGA_sur_single_file == "pdf" & input$TCGA_sur_single_width > 50){
        w = 50
      }else if(input$TCGA_sur_single_file != "pdf" & input$TCGA_sur_single_width <= 3000){
        w = input$TCGA_sur_single_width
      }else{
        w = 3000
      }
      
      if(input$TCGA_sur_single_file == 'pdf'){
        pdf(file = file,width = w,height = h)
      }else if(input$TCGA_sur_single_file == 'png'){
        png(file = file,width = w,height = h,res = r)
      }else if(input$TCGA_sur_single_file == 'jpeg'){
        jpeg(file = file,width = w,height = h,res = r)
      }else{
        tiff(file = file,width = w,height = h,res = r)
      }
      
      
      print(TCGA_survival(TCGA = TCGA,cancer_type = cancer_type,gene = gene,mut_patient_id = TCGA_cohort$mut,wt_patient_id = TCGA_cohort$wt))
      
      dev.off()
    }
  )
  
}