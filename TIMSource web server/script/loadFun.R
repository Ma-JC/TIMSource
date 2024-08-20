library("shinythemes")
library("DT")
library("survival")
library("survminer")
library("scales")
library("shinydashboard")
library("maftools")
library("patchwork")
library("plotly")
library("shinycssloaders")
library("limma")
library("reshape2")
library("shinyvalidate")
library("shinyWidgets")
library("gridExtra")
library("enrichplot")
library("clusterProfiler")
library("RColorBrewer")
library("shinyBS")
library("cicerone")
library("grid")
library("reactablefmtr")

library("forestplot")
library("ComplexHeatmap")
library("httr")
library("rvest")
library("bsplus")
library("callr")
library("waiter")
library("installr")
library("meta")
myGSEA = function(geneList,TERM2GENE,pvalueCutoff){
  library(clusterProfiler)
  GSEA(geneList = geneList,TERM2GENE = TERM2GENE,pvalueCutoff = pvalueCutoff)
}
##################Ref_datasets##########################

ref_cohort_cal = function(dataset,gene,dataset_mu,Mut_type,Wild_type,therapy_type){
  
  tp = rownames(dataset)[dataset$THERAPY %in% therapy_type]
  
    mut = unique(dataset_mu$ID[
      dataset_mu$Hugo_Symbol == gene & 
        dataset_mu$Variant_Classification %in% Mut_type
    ])
  
  if(Wild_type == "Others"){
    wt = setdiff(unique(c(dataset_mu$ID,rownames(dataset))),mut)
  }else{
    wt = setdiff(unique(c(dataset_mu$ID,rownames(dataset))),unique(dataset_mu$ID[dataset_mu$Hugo_Symbol == gene]))
  }
  
    mut = intersect(tp,mut)
    wt = intersect(tp,wt)
  return(list("mut"=mut,"wt"=wt))
  
}

os_survival = function(dataset,gene,dataset_mu,mut,wt){
  
  tmp_data = dataset[,c("OS_TIME","OS_STATUS")]
  
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  
  fit <- do.call(survfit, list(Surv(OS_TIME,OS_STATUS)~groups,data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "Month",
                  pval = T,
                  pval.size = 6,
                  fontsize = 6,
                  font.x = c(18, 'plain', 'black'),
                  font.y = c(18, 'plain', 'black'), 
                  risk.table = F,
                  cumevents.y.text = FALSE, 
                  cumevents.y.text.col = TRUE,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival(OS)",
                  legend.title=gene,
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
                                                    legend.title = element_text(size = 16,face = "bold"),
                                                    legend.text = element_text(size = 16))
  )
  return(p1)
  
  
  
}

pfs_survival = function(dataset,gene,dataset_mu,mut,wt){
  
  tmp_data = dataset[,c("PFS_TIME","PFS_STATUS")]
  
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  
  fit <- do.call(survfit, list(Surv(PFS_TIME,PFS_STATUS)~groups,data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "Month",
                  pval = T,
                  pval.size = 6,
                  fontsize = 6,
                  font.x = c(18, 'plain', 'black'), 
                  font.y = c(18, 'plain', 'black'), 
                  risk.table = F,
                  cumevents.y.text = FALSE, 
                  cumevents.y.text.col = TRUE,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Progress Free Survival(PFS)",
                  legend.title=gene,
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
                                                    legend.title = element_text(size = 16,face = "bold"),
                                                    legend.text = element_text(size = 16))
  )
  return(p1)
  
  
}

CRPR = function(dataset,gene,dataset_mu,mut,wt){
  
    tmp_data = dataset["RECIST"]
    
    tmp_data$groups = NA
    tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
    tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
    
    contingency = table(tmp_data$groups,tmp_data$RECIST)
    rc = rowSums(contingency)
    tmp_ma = as.data.frame.matrix(contingency/rc)
    

  

	if(sum(contingency)<40){
	  pvalue = fisher.test(contingency)$p.value
	}else{
	  tmp_chi = chisq.test(contingency,correct = F)
	  if(sum(tmp_chi$expected < 1) > 0){
		pvalue = fisher.test(contingency)$p.value
	  }else if(sum(tmp_chi$expected < 5) > 0){
		pvalue = chisq.test(contingency,correct = TRUE)$p.value
	  }else{pvalue = chisq.test(contingency,correct = FALSE)$p.value}
	}


	values = c(tmp_ma[,1],tmp_ma[,2])
	g1 = c("Mutation","Wildtype","Mutation","Wildtype")
	g2 = c("CR/PR","CR/PR","SD/PD","SD/PD")
	data = data.frame("Genetype" = g1,"CR_PR_ratio" = values,"response" = g2)
	data$response = factor(data$response,levels = c("SD/PD","CR/PR"))

	p = ggplot(data = data)+
  	  ggtitle("Drugs Response (RECIST)")+
  	  geom_bar(mapping = aes(x = Genetype,y = CR_PR_ratio,fill = response),stat = 'identity', position = 'stack')+
  	  geom_text(aes(x = Genetype,y = CR_PR_ratio,group = response,label = round(CR_PR_ratio,4)),size = 10,stat = "identity", position = position_stack(0.5))+
  	  theme_classic()+
  	  scale_y_continuous(expand = c(0,0),limits = c(0,1.15),breaks = c(0,0.25,0.5,0.75,1))+
  	  scale_fill_manual(values = c("sandybrown","skyblue1"))+
  	  annotate("text",x = 1.5,y = 1.1,label = paste("Wilcox : P value = ",round(pvalue,digits = 4)),size = 6)+
  	  labs(fill = gene)+
  	  ylab("Patient")+
  	  theme(
  		plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
  		axis.text.x = element_text(size = 15,face = "bold"),
  		axis.title.x.bottom = element_blank(),
  		axis.title.y = element_text(size = 15,face = "bold"),
  		axis.text.y.left = element_text(size = 15),
  		legend.text = element_text(size = 15),
  		legend.title = element_text(size = 15,face = "bold"),
  		legend.position = "right")
	
	return(p)
}

DCB = function(dataset,gene,dataset_mu,mut,wt){
  
    tmp_data = dataset["RESPONSE"]
    
    tmp_data$groups = NA
    tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
    tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
    
    contingency = table(tmp_data$groups,tmp_data$RESPONSE)
    rc = rowSums(contingency)
    tmp_ma = as.data.frame.matrix(contingency/rc)
    
  
  


	if(sum(contingency)<40){
	  pvalue = fisher.test(contingency)$p.value
	}else{
	  tmp_chi = chisq.test(contingency,correct = F)
	  if(sum(tmp_chi$expected < 1) > 0){
		pvalue = fisher.test(contingency)$p.value
	  }else if(sum(tmp_chi$expected < 5) > 0){
		pvalue = chisq.test(contingency,correct = TRUE)$p.value
	  }else{pvalue = chisq.test(contingency,correct = FALSE)$p.value}
	}


	values = c(tmp_ma[,1],tmp_ma[,2])
	g1 = c("Mutation","Wildtype","Mutation","Wildtype")
	g2 = c("NCB","NCB","DCB","DCB")
	data = data.frame("Genetype" = g1,"CR_PR_ratio" = values,"response" = g2)
	data$response = factor(data$response,levels = c("NCB","DCB"))

	p = ggplot(data = data)+
  	  ggtitle("Clinical Benefit")+
  	  geom_bar(mapping = aes(x = Genetype,y = CR_PR_ratio,fill = response),stat = 'identity', position = 'stack')+
  	  geom_text(aes(x = Genetype,y = CR_PR_ratio,group = response,label = round(CR_PR_ratio,4)),size = 10,stat = "identity", position = position_stack(0.5))+
  	  theme_classic()+
  	  scale_y_continuous(expand = c(0,0),limits = c(0,1.15),breaks = c(0,0.25,0.5,0.75,1))+
  	  scale_fill_manual(values = c("sandybrown","skyblue1"))+
  	  annotate("text",x = 1.5,y = 1.1,label = paste("Wilcox : P value = ",round(pvalue,digits = 4)),size = 6)+
  	  labs(fill = gene)+
  	  ylab("Patient")+
  	  theme(
  		plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
  		axis.text.x = element_text(size = 15,face = "bold"),
  		axis.title.x.bottom = element_blank(),
  		axis.title.y = element_text(size = 15,face = "bold"),
  		axis.text.y.left = element_text(size = 15),
  		legend.text = element_text(size = 15),
  		legend.title = element_text(size = 15,face = "bold"),
  		legend.position = "right")
	
	return(p)
}

TMB = function(dataset,gene,remove_outlier = TRUE,dataset_mu,mut,wt){

  tmp_data = dataset["TMB"]

  
  tmp_data$gene = NA
  tmp_data$gene[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$gene[rownames(tmp_data) %in% wt] = "Wildtype"
  
  tmp_data = tmp_data[!is.na(tmp_data$gene),]
  res = wilcox.test(tmp_data[,"TMB"]~tmp_data[,"gene"],data = tmp_data)
  if(remove_outlier){
    lim1 = boxplot.stats(tmp_data[tmp_data$gene == "Mutation","TMB"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$gene == "Wildtype","TMB"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = gene,y = TMB,color = gene),alpha = 0.5,outlier.color = NA)+
      ggtitle("Tumor Mutation Burden (TMB)")+
      labs(color=gene)+
      geom_point(mapping = aes(x = gene,y = TMB,color = gene),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "top"
      ) +
      scale_color_manual(values = c("red2", "blue4"))
    
    return(p)
  }else{
    max_val = max(tmp_data$TMB)
    min_val = min(tmp_data$TMB)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = gene,y = TMB,color = gene),alpha = 0.5,outlier.colour = NA)+
        ggtitle("Tumor Mutation Burden (TMB)")+
        labs(color=gene)+
        geom_point(mapping = aes(x = gene,y = TMB,color = gene),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        scale_color_manual(values = c("red2", "blue4"))+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    
    return(p)
    
  }
  
}

MutationType = function(dataset_mu,gene){
  tmp_data = table(dataset_mu[dataset_mu$Hugo_Symbol %in% gene,"Variant_Classification"])
  tmp_data = tmp_data/sum(tmp_data)
  tmp = data.frame("Mutation_Type" = names(tmp_data),"values" = as.numeric(tmp_data))
  p = ggplot(data = tmp)+
      ggtitle("Variant_Classification")+
      geom_bar(mapping = aes(x = Mutation_Type,y = values,fill = Mutation_Type),stat = "identity")+
      theme_classic()+
      coord_flip()+
      scale_y_continuous(expand = c(0,0),limits = c(0,1))+
      theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "right")   
  return(p)
  
}

DEG_ref = function(datasets_rna_wes,dataset ,min.pct,dataset_mu,mut,wt){
  
  
  mut = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),mut)
  wt = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),wt)
  
  tmp_data = datasets_rna_wes[[dataset]][["RNA"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

immune_one_datasets = function(datasets_rna_wes,dataset,gene,immune_module,selection,outlier,dataset_mu,mut,wt){

  
  mut = intersect(colnames(datasets_rna_wes[[dataset]][[immune_module]]),mut)
  wt = intersect(colnames(datasets_rna_wes[[dataset]][[immune_module]]),wt)
  
  tmp_data = data.frame(values= datasets_rna_wes[[dataset]][[immune_module]][selection,c(mut,wt)],groups = c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))))
  res = wilcox.test(values~groups,data = tmp_data)
  
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "Mutation","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Wildtype","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(selection)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
              plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
              axis.text.x = element_text(size = 15,face = "bold"),
              axis.title.x.bottom = element_blank(),
              axis.title.y = element_text(size = 15,face = "bold"),
              axis.text.y.left = element_text(size = 15),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15,face = "bold"),
              legend.position = "top")
    return(p)
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
        ggtitle(selection)+
        labs(color=gene)+
        ylab("Score")+
        geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        scale_color_manual(values = c("red2", "blue4"))+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    return(p)
    
  } 
  
}

immune_all_datasets = function(datasets_rna_wes,dataset,gene,immune_module,dataset_mu,mut,wt){
  
  
  mut = intersect(colnames(datasets_rna_wes[[dataset]][[immune_module]]),mut)
  wt = intersect(colnames(datasets_rna_wes[[dataset]][[immune_module]]),wt)
  
  wt_values = datasets_rna_wes[[dataset]][[immune_module]][,wt]
  mut_values = datasets_rna_wes[[dataset]][[immune_module]][,mut]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune"
  wt_values$groups = "Wiltype"
  
  if(immune_module == "immune_infiltration"){
    immune_title = "Infiltrating immune cells"
  }else{immune_title = "Immune-related signatures"}
  tmp_data = rbind(mut_values,wt_values)
  # tmp_data$immune = as.factor(tmp_data$immune)
  p = ggplot(data = tmp_data)+
    ggtitle(immune_title)+
    geom_boxplot(aes(x = immune,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune == i,])
    max_val = max(tmp_data[tmp_data$immune == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

##################Ref_datasets_pm##########################

ref_cohort_cal_pm = function(dataset,gene,dataset_mu,Mut_type,Wild_type,therapy_type){
  
  tp = rownames(dataset)[dataset$THERAPY %in% therapy_type]
  
  genes = intersect(colnames(dataset),pathway_list[[gene]])

    
    
    mut = unique(dataset_mu$ID[
      dataset_mu$Hugo_Symbol %in% genes & 
        dataset_mu$Variant_Classification %in% Mut_type
    ])
  
  if(Wild_type == "Others"){
    wt = setdiff(unique(c(dataset_mu$ID,rownames(dataset))),mut)
  }else{
    wt = setdiff(unique(c(dataset_mu$ID,rownames(dataset))),unique(dataset_mu$ID[dataset_mu$Hugo_Symbol %in% genes]))
  }
  
  mut = intersect(tp,mut)
  wt = intersect(tp,wt)
  return(list("mut"=mut,"wt"=wt))
  
}

os_survival_pm = function(dataset,gene,dataset_mu,mut,wt){
  

  tmp_data = dataset[,c("OS_TIME","OS_STATUS")]
  
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  
  fit <- do.call(survfit, list(Surv(OS_TIME,OS_STATUS)~groups,data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "Month",
                  pval = T,
                  pval.size = 6,
                  fontsize = 6,
                  font.x = c(18, 'plain', 'black'),
                  font.y = c(18, 'plain', 'black'), 
                  risk.table = F,
                  cumevents.y.text = FALSE, 
                  cumevents.y.text.col = TRUE,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival(OS)",
                  legend.title="Pathway",
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
                                                    legend.title = element_text(size = 16,face = "bold"),
                                                    legend.text = element_text(size = 16))
  )
  tmp_table = p1$cumevents + theme(axis.title.y = element_blank())
  p1$cumevents = tmp_table
  return(p1)
  
  
  

}

pfs_survival_pm = function(dataset,gene,dataset_mu,mut,wt){

  tmp_data = dataset[,c("PFS_TIME","PFS_STATUS")]
  
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  
  fit <- do.call(survfit, list(Surv(PFS_TIME,PFS_STATUS)~groups,data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "Month",
                  pval = T,
                  pval.size = 6,
                  fontsize = 6,
                  font.x = c(18, 'plain', 'black'), 
                  font.y = c(18, 'plain', 'black'), 
                  risk.table = F,
                  cumevents.y.text = FALSE, 
                  cumevents.y.text.col = TRUE,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Progress Free Survival(PFS)",
                  legend.title="Pathway",
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
                                                    legend.title = element_text(size = 16,face = "bold"),
                                                    legend.text = element_text(size = 16))
  )
  tmp_table = p1$cumevents + theme(axis.title.y = element_blank())
  p1$cumevents = tmp_table
  return(p1)
}

CRPR_pm = function(dataset,gene,dataset_mu,mut,wt){
  
  tmp_data = dataset[c("RECIST")]
  
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  

  contingency = table(tmp_data$groups,tmp_data$RECIST)
  rc = rowSums(contingency)
  tmp_ma = as.data.frame.matrix(contingency/rc)
  
  if(sum(contingency)<40){
    pvalue = fisher.test(contingency)$p.value
  }else{
    tmp_chi = chisq.test(contingency,correct = F)
    if(sum(tmp_chi$expected < 1) > 0){
      pvalue = fisher.test(contingency)$p.value
    }else if(sum(tmp_chi$expected < 5) > 0){
      pvalue = chisq.test(contingency,correct = TRUE)$p.value
    }else{pvalue = chisq.test(contingency,correct = FALSE)$p.value}
  }
  
  
  values = c(tmp_ma[,1],tmp_ma[,2])
  g1 = c("Mutation","Wildtype","Mutation","Wildtype")
  g2 = c("CR/PR","CR/PR","SD/PD","SD/PD")
  data = data.frame("Genetype" = g1,"CR_PR_ratio" = values,"response" = g2)
  data$response = factor(data$response,levels = c("SD/PD","CR/PR"))
  
  p = ggplot(data = data)+
      ggtitle("Drugs Response (RECIST)")+
      geom_bar(mapping = aes(x = Genetype,y = CR_PR_ratio,fill = response),stat = 'identity', position = 'stack')+
      geom_text(aes(x = Genetype,y = CR_PR_ratio,group = response,label = round(CR_PR_ratio,4)),size = 10,stat = "identity", position = position_stack(0.5))+
      theme_classic()+
      scale_y_continuous(expand = c(0,0),limits = c(0,1.15),breaks = c(0,0.25,0.5,0.75,1))+
      scale_fill_manual(values = c("sandybrown","skyblue1"))+
      annotate("text",x = 1.5,y = 1.1,label = paste("Wilcox : P value = ",round(pvalue,digits = 4)),size = 6)+
      labs(fill = "PATHWAY")+
      ylab("Patient")+
      theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "right")
  
  return(p)
}

DCB_pm = function(dataset,gene,dataset_mu,mut,wt){
  

  tmp_data = dataset[c("RESPONSE")]
  
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  

  contingency = table(tmp_data$groups,tmp_data$RESPONSE)
  rc = rowSums(contingency)
  tmp_ma = as.data.frame.matrix(contingency/rc)
  
  if(sum(contingency)<40){
    pvalue = fisher.test(contingency)$p.value
  }else{
    tmp_chi = chisq.test(contingency,correct = F)
    if(sum(tmp_chi$expected < 1) > 0){
      pvalue = fisher.test(contingency)$p.value
    }else if(sum(tmp_chi$expected < 5) > 0){
      pvalue = chisq.test(contingency,correct = TRUE)$p.value
    }else{pvalue = chisq.test(contingency,correct = FALSE)$p.value}
  }
  
  
  values = c(tmp_ma[,1],tmp_ma[,2])
  g1 = c("Mutation","Wildtype","Mutation","Wildtype")
  g2 = c("NCB","NCB","DCB","DCB")
  data = data.frame("Genetype" = g1,"CR_PR_ratio" = values,"response" = g2)
  data$response = factor(data$response,levels = c("NCB","DCB"))
  
  p = ggplot(data = data)+
      ggtitle("Clinical Benefit")+
      geom_bar(mapping = aes(x = Genetype,y = CR_PR_ratio,fill = response),stat = 'identity', position = 'stack')+
      geom_text(aes(x = Genetype,y = CR_PR_ratio,group = response,label = round(CR_PR_ratio,4)),size = 10,stat = "identity", position = position_stack(0.5))+
      theme_classic()+
      scale_y_continuous(expand = c(0,0),limits = c(0,1.15),breaks = c(0,0.25,0.5,0.75,1))+
      scale_fill_manual(values = c("sandybrown","skyblue1"))+
      annotate("text",x = 1.5,y = 1.1,label = paste("Wilcox : P value = ",round(pvalue,digits = 4)),size = 6)+
      labs(fill = "PATHWAY")+
      ylab("Patient")+
      theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "right")
  
  return(p)
}

TMB_pm = function(dataset,gene,remove_outlier = TRUE,dataset_mu,mut,wt){
  

  tmp_data = dataset[c("TMB")]
  
  tmp_data$gene = NA
  tmp_data$gene[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$gene[rownames(tmp_data) %in% wt] = "Wildtype"
  tmp_data = tmp_data[!is.na(tmp_data$gene),]

  res = wilcox.test(tmp_data[,"TMB"]~tmp_data[,"gene"],data = tmp_data)
  if(remove_outlier){
    lim1 = boxplot.stats(tmp_data[tmp_data$gene == "Mutation","TMB"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$gene == "Wildtype","TMB"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = gene,y = TMB,color = gene),alpha = 0.5,outlier.color = NA)+
      ggtitle("Tumor Mutation Burden (TMB)")+
      labs(color="PATHWAY")+
      geom_point(mapping = aes(x = gene,y = TMB,color = gene),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "top"
        ) +
      scale_color_manual(values = c("red2", "blue4"))
    return(p)
    
  }else{
    max_val = max(tmp_data$TMB)
    min_val = min(tmp_data$TMB)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = gene,y = TMB,color = gene),alpha = 0.5,outlier.colour = NA)+
        ggtitle("Tumor Mutation Burden (TMB)")+
        labs(color=gene)+
        geom_point(mapping = aes(x = gene,y = TMB,color = gene),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")+
      scale_color_manual(values = c("red2", "blue4"))
    
    return(p)
  }
  
}

MutationType_pm = function(dataset_mu,gene){
  gene = intersect(dataset_mu$Hugo_Symbol,pathway_list[[gene]])
  tmp_data = table(dataset_mu[dataset_mu$Hugo_Symbol %in% gene,"Variant_Classification"])
  tmp_data = tmp_data/sum(tmp_data)
  tmp = data.frame("Mutation_Type" = names(tmp_data),"values" = as.numeric(tmp_data))
  p = ggplot(data = tmp)+
      ggtitle("Variant_Classification")+
      geom_bar(mapping = aes(x = Mutation_Type,y = values,fill = Mutation_Type),stat = "identity")+
      theme_classic()+
      coord_flip()+
      scale_y_continuous(expand = c(0,0),limits = c(0,1))+
      theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "right")   
  return(p)
  
  
}

DEG_ref_pm = function(datasets_rna_wes,dataset,gene ,min.pct,dataset_mu,mut,wt){
  
  mut = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),mut)
  wt = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),wt)
  
  tmp_data = datasets_rna_wes[[dataset]][["RNA"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

immune_one_datasets_pm = function(datasets_rna_wes,dataset,gene,immune_module,selection,outlier,dataset_mu,mut,wt){
  
  mut = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),mut)
  wt = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),wt)
  

  tmp_data = data.frame(values= datasets_rna_wes[[dataset]][[immune_module]][selection,c(mut,wt)],groups = c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))))
  res = wilcox.test(values~groups,data = tmp_data)
  
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "Mutation","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Wildtype","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(selection)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
            plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
            axis.text.x = element_text(size = 15,face = "bold"),
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 15,face = "bold"),
            axis.text.y.left = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15,face = "bold"),
            legend.position = "top")
    
    return(p)
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
        ggtitle(selection)+
        labs(color=gene)+
        ylab("Score")+
        scale_color_manual(values = c("red2", "blue4"))+
        geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    
    return(p)
    
  } 
  
}

immune_all_datasets_pm = function(datasets_rna_wes,dataset,gene,immune_module,dataset_mu,mut,wt){
  
  mut = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),mut)
  wt = intersect(colnames(datasets_rna_wes[[dataset]]$RNA),wt)
  
  wt_values = datasets_rna_wes[[dataset]][[immune_module]][,wt]
  mut_values = datasets_rna_wes[[dataset]][[immune_module]][,mut]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune"
  wt_values$groups = "Wiltype"
  
  if(immune_module == "immune_infiltration"){
    immune_title = "Infiltrating immune cells"
  }else{immune_title = "Immune-related signatures"}
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle(immune_title)+
    geom_boxplot(aes(x = immune,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = "PATHWAY")+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune == i,])
    max_val = max(tmp_data[tmp_data$immune == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}
#######################TCGA################################
TCGA_cohort_cal = function(TCGA,cancer_type,gene,Mut_type,Wild_type){
  if("All" %in% Mut_type){
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol == gene
    ]))
  }else{
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol == gene &
        TCGA[[cancer_type]]$maf@data$Variant_Classification %in% Mut_type
    ]))
  }
  
  if(Wild_type == "Others"){
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode)),mut) 
  }else{
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode)),
                 unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[TCGA[[cancer_type]]$maf@data$Hugo_Symbol == gene])))
  }
  return(list("mut" = mut,"wt" = wt))
}

Mutational_Landscape = function(TCGA,cancer_type,gene,mut,wt){
  
  tmp1 <- subsetMaf(maf=TCGA[[cancer_type]][["maf"]], tsb=mut, isTCGA=TRUE)
  tmp2 <- subsetMaf(maf=TCGA[[cancer_type]][["maf"]], tsb=wt, isTCGA=TRUE)
  
  fvsm <- mafCompare(m1=tmp1, m2=tmp2, m1Name="Mutation", m2Name="Wildtype", minMut=5)
  
  meta = data.frame("Tumor_Sample_Barcode" = c(mut,wt),"anno" = c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))))
  colnames(meta)[2] = gene
  
  return(list("fvsm"=fvsm,"meta"=meta))
}

immune_one = function(TCGA,cancer_type,gene,immune_module,selection,outlier,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
  wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
  
  tmp_data = data.frame(values=c(mut_values,wt_values),groups = c(rep("Mutation",length(mut_values)),rep("Wildtype",length(wt_values))))
  res = wilcox.test(values~groups,data = tmp_data)
  
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "Mutation","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Wildtype","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(selection)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p1+      theme(
      plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
      axis.text.x = element_text(size = 15,face = "bold"),
      axis.title.x.bottom = element_blank(),
      axis.title.y = element_text(size = 15,face = "bold"),
      axis.text.y.left = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15,face = "bold"),
      legend.position = "top")
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
      ggtitle(selection)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
      theme_classic()+
      annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
      theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "top")
    
  } 
  
  
}

immune_infiltration = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_infiltration"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_infiltration"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("infiltrating immune cells")+
    geom_boxplot(aes(x = immune_infiltration,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_pathway"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_pathway"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures")+
    geom_boxplot(aes(x = immune_pathway,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

DEG = function(TCGA,gene,cancer_type,min.pct,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["rna"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["rna"]]),wt_patient_id)
  
  tmp_data = TCGA[[cancer_type]][["rna"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

TCGA_survival = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),mut_patient_id)
  wt = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),wt_patient_id)
  tmp_data = TCGA[[cancer_type]][["clinical_detial"]]
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  tmp_data = tmp_data[,-1]
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title=gene,
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  
  
  if(cancer_type != "LAML"){
    fit <- do.call(survfit, list(Surv(DSS.time,DSS)~tmp_data[,"groups"],data = tmp_data))
    p2 = ggsurvplot(fit,data = tmp_data,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Disease-specific survival (DSS)",
                    legend.title=gene,
                    legend.labs = c("Mutation","Wildtype"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    
    fit <- do.call(survfit, list(Surv(PFI.time,PFI)~tmp_data[,"groups"],data = tmp_data))
    p3 = ggsurvplot(fit,data = tmp_data,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Progression-free interval (PFI)",
                    legend.title=gene,
                    legend.labs = c("Mutation","Wildtype"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))

  }
  
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  
  if(cancer_type != "LAML"){
    p = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)  | (p3$plot/p3$cumevents)+plot_layout(design = layout)
    return(p)
  }else{
    p = (p1$plot/p1$cumevents)+plot_layout(design = layout) 
    return(p)
  }
  
  
}

#######################TCGA_pm################################
TCGA_cohort_cal_pm = function(TCGA,cancer_type,gene,Mut_type,Wild_type){
  if("All" %in% Mut_type){
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol %in% pathway_list[[gene]]
    ]))
  }else{
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol %in% pathway_list[[gene]] &
        TCGA[[cancer_type]]$maf@data$Variant_Classification %in% Mut_type
    ]))
  }
  
  if(Wild_type == "Others"){
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode)),mut) 
  }else{
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode)),
                 unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[TCGA[[cancer_type]]$maf@data$Hugo_Symbol %in% pathway_list[[gene]]])))
  }
  return(list("mut" = mut,"wt" = wt))
}

Mutational_Landscape_pm = function(TCGA,cancer_type,gene,mut,wt){
  
  tmp1 <- subsetMaf(maf=TCGA[[cancer_type]][["maf"]], tsb=mut, isTCGA=TRUE)
  tmp2 <- subsetMaf(maf=TCGA[[cancer_type]][["maf"]], tsb=wt, isTCGA=TRUE)
  
  fvsm <- mafCompare(m1=tmp1, m2=tmp2, m1Name="Mutation", m2Name="Wildtype", minMut=5)
  
  meta = data.frame("Tumor_Sample_Barcode" = c(mut,wt),"anno" = c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))))
  colnames(meta)[2] = gene
  
  tmp_maf = subsetMaf(maf = TCGA[[cancer_type]][["maf"]],tsb = meta$Tumor_Sample_Barcode)
  
  genes = TCGA[[cancer_type]]$maf@gene.summary$Hugo_Symbol[TCGA[[cancer_type]]$maf@gene.summary$Hugo_Symbol %in% pathway_list[[gene]]]
  if(length(genes) > 50){
    genes = genes[1:50]
  }
  return(list("fvsm"=fvsm,"meta"=meta,"tmp_maf" = tmp_maf,"genes" = genes))
}

immune_one_pm = function(TCGA,cancer_type,gene,immune_module,selection,outlier,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
  wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
  
  tmp_data = data.frame(values=c(mut_values,wt_values),groups = c(rep("Mutation",length(mut_values)),rep("Wildtype",length(wt_values))))
  res = wilcox.test(values~groups,data = tmp_data)
  
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "Mutation","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Wildtype","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(selection)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p1+      theme(
      plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
      axis.text.x = element_text(size = 15,face = "bold"),
      axis.title.x.bottom = element_blank(),
      axis.title.y = element_text(size = 15,face = "bold"),
      axis.text.y.left = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15,face = "bold"),
      legend.position = "top")
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
      ggtitle(selection)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
      theme_classic()+
      annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
      theme(
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold"),
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.y.left = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15,face = "bold"),
        legend.position = "top")
    
  } 
  
  
}

immune_infiltration_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_infiltration"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_infiltration"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("infiltrating immune cells")+
    geom_boxplot(aes(x = immune_infiltration,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_pathway"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_pathway"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures")+
    geom_boxplot(aes(x = immune_pathway,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

DEG_pm = function(TCGA,gene,cancer_type,min.pct,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["rna"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["rna"]]),wt_patient_id)
  
  tmp_data = TCGA[[cancer_type]][["rna"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

TCGA_survival_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),mut_patient_id)
  wt = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),wt_patient_id)
  tmp_data = TCGA[[cancer_type]][["clinical_detial"]]
  tmp_data$groups = NA
  tmp_data$groups[rownames(tmp_data) %in% mut] = "Mutation"
  tmp_data$groups[rownames(tmp_data) %in% wt] = "Wildtype"
  tmp_data = tmp_data[,-1]
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title=gene,
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  
  
  
  if(cancer_type != "LAML"){
    fit <- do.call(survfit, list(Surv(DSS.time,DSS)~tmp_data[,"groups"],data = tmp_data))
    p2 = ggsurvplot(fit,data = tmp_data,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Disease-specific survival (DSS)",
                    legend.title=gene,
                    legend.labs = c("Mutation","Wildtype"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    
    fit <- do.call(survfit, list(Surv(PFI.time,PFI)~tmp_data[,"groups"],data = tmp_data))
    p3 = ggsurvplot(fit,data = tmp_data,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Progression-free interval (PFI)",
                    legend.title=gene,
                    legend.labs = c("Mutation","Wildtype"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  }

  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  if(cancer_type != "LAML"){
    p = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)  | (p3$plot/p3$cumevents)+plot_layout(design = layout)
    return(p)
  }else{
    p = (p1$plot/p1$cumevents)+plot_layout(design = layout) 
    return(p)
  }  
  
}

#######################CPTAC################################
CPTAC_cohort_cal = function(TCGA,cancer_type,gene,Mut_type,Wild_type){
  if("All" %in% Mut_type){
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol == gene
    ]))
  }else{
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol == gene &
        TCGA[[cancer_type]]$maf@data$Variant_Classification %in% Mut_type
    ]))
  }
  
  if(Wild_type == "Others"){
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA)),mut) 
  }else{
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA)),
                 unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[TCGA[[cancer_type]]$maf@data$Hugo_Symbol == gene])))
  }
  return(list("mut" = mut,"wt" = wt))
}

Mutational_Landscape_cptac = function(TCGA,cancer_type,gene,mut,wt){
  
  tmp = TCGA[[cancer_type]][["maf"]]
  tmp@data$Tumor_Sample_Barcode = tmp@data$TCGA
  tmp@data$Tumor_Sample_Barcode = as.factor(tmp@data$Tumor_Sample_Barcode)
  tmp1 <- subsetMaf(maf=tmp,  tsb=mut, isTCGA=T)
  tmp2 <- subsetMaf(maf=tmp,  tsb=wt, isTCGA=T)
  
  tmp1@summary$summary[3] = length(mut)
  tmp2@summary$summary[3] = length(wt)
  
  fvsm <- mafCompare(m1=tmp1, m2=tmp2, m1Name="Mutation", m2Name="Wildtype", minMut=5)
  
  meta = data.frame("Tumor_Sample_Barcode" = c(mut,wt),"anno" = c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))))
  colnames(meta)[2] = gene
  
  tmp_maf = subsetMaf(maf = tmp,tsb = meta$Tumor_Sample_Barcode)
  tmp_maf@summary$summary[3] = length(wt) + length(mut)
  
  return(list("fvsm"=fvsm,"meta"=meta,"tmp_maf" = tmp_maf))
}

immune_one_CPTAC = function(TCGA,cancer_type,gene,immune_module,selection,outlier,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
  wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
  
  tmp_data = data.frame(values=c(mut_values,wt_values),groups = c(rep("Mutation",length(mut_values)),rep("Wildtype",length(wt_values))))
  res = wilcox.test(values~groups,data = tmp_data)
  if(grepl("protein",immune_module)){title1 = paste(selection,"(Protein)")}else{title1 = paste(selection,"(RNA)")}
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "Mutation","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Wildtype","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(title1)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
            plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
            axis.text.x = element_text(size = 15,face = "bold"),
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 15,face = "bold"),
            axis.text.y.left = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15,face = "bold"),
            legend.position = "top")
    return(p)
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
        ggtitle(title1)+
        labs(color=gene)+
        ylab("Score")+
        scale_color_manual(values = c("red2", "blue4"))+
        geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    return(p)
  } 
  
  
}

immune_infiltration_CPTAC_rna = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_infiltration"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_infiltration"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(RNA)")+
    geom_boxplot(aes(x = immune_infiltration,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_CPTAC_protein = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_infiltration_protein"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_infiltration_protein"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(protein)")+
    geom_boxplot(aes(x = immune_infiltration_protein,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration_protein)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration_protein == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration_protein == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_CPTAC_rna = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_pathway"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_pathway"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (RNA)")+
    geom_boxplot(aes(x = immune_pathway,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_CPTAC_protein = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_pathway_protein"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_pathway_protein"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (protein)")+
    geom_boxplot(aes(x = immune_pathway_protein,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway_protein)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway_protein == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway_protein == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

DEG_CPTAC_rna = function(TCGA,gene,cancer_type,min.pct,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["rna"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["rna"]]),wt_patient_id)
  
  tmp_data = TCGA[[cancer_type]][["rna"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

DEG_CPTAC_protein = function(TCGA,gene,cancer_type,min.pct,mut_patient_id,wt_patient_id){

  mut = intersect(colnames(TCGA[[cancer_type]][["protein"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["protein"]]),wt_patient_id)
  
  tmp_data = TCGA[[cancer_type]][["protein"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

CPTAC_survival = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,mut_patient_id)
  wt = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,wt_patient_id)
  tmp_data = TCGA[[cancer_type]][["Survival"]]
  tmp_data$groups = NA
  tmp_data$groups[tmp_data$PatientID %in% mut] = "Mutation"
  tmp_data$groups[tmp_data$PatientID %in% wt] = "Wildtype"
  tmp_data = tmp_data[,-1]
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title=gene,
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  p = (p1$plot/p1$cumevents)+plot_layout(design = layout)
  return(p)
  
}

#######################CPTAC_pm################################
CPTAC_cohort_cal_pm = function(TCGA,cancer_type,gene,Mut_type,Wild_type){
  if("All" %in% Mut_type){
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol %in% pathway_list[[gene]]
    ]))
  }else{
    mut = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[
      TCGA[[cancer_type]]$maf@data$Hugo_Symbol %in% pathway_list[[gene]] &
        TCGA[[cancer_type]]$maf@data$Variant_Classification %in% Mut_type
    ]))
  }
  
  if(Wild_type == "Others"){
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA)),mut) 
  }else{
    wt = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA)),
                 unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[TCGA[[cancer_type]]$maf@data$Hugo_Symbol %in% pathway_list[[gene]]])))
  }
  return(list("mut" = mut,"wt" = wt))
}

Mutational_Landscape_cptac_pm = function(TCGA,cancer_type,gene,mut,wt){

  
  tmp = TCGA[[cancer_type]][["maf"]]
  tmp@data$Tumor_Sample_Barcode = tmp@data$TCGA
  tmp@data$Tumor_Sample_Barcode = as.factor(tmp@data$Tumor_Sample_Barcode)
  tmp1 <- subsetMaf(maf=tmp,  tsb=mut, isTCGA=T)
  tmp2 <- subsetMaf(maf=tmp,  tsb=wt, isTCGA=T)
  
  tmp1@summary$summary[3] = length(mut)
  tmp2@summary$summary[3] = length(wt)
  
  fvsm <- mafCompare(m1=tmp1, m2=tmp2, m1Name="Mutation", m2Name="Wildtype", minMut=5)
  
  meta = data.frame("Tumor_Sample_Barcode" = c(mut,wt),"anno" = c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))))
  colnames(meta)[2] = gene
  
  tmp_maf = subsetMaf(maf = tmp,tsb = meta$Tumor_Sample_Barcode)
  tmp_maf@summary$summary[3] = length(wt) + length(mut)

  genes = TCGA[[cancer_type]]$maf@gene.summary$Hugo_Symbol[TCGA[[cancer_type]]$maf@gene.summary$Hugo_Symbol %in% pathway_list[[gene]]]
  if(length(genes) > 50){
    genes = genes[1:50]
  }
  return(list("fvsm"=fvsm,"meta"=meta,"tmp_maf" = tmp_maf,"genes" = genes))
}

immune_one_CPTAC_pm = function(TCGA,cancer_type,gene,immune_module,selection,outlier,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
  wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
  
  tmp_data = data.frame(values=c(mut_values,wt_values),groups = c(rep("Mutation",length(mut_values)),rep("Wildtype",length(wt_values))))
  res = wilcox.test(values~groups,data = tmp_data)
  if(grepl("protein",immune_module)){title1 = paste(selection,"(Protein)")}else{title1 = paste(selection,"(RNA)")}
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "Mutation","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Wildtype","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(title1)+
      labs(color=gene)+
      ylab("Score")+
      scale_color_manual(values = c("red2", "blue4"))+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
            plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
            axis.text.x = element_text(size = 15,face = "bold"),
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 15,face = "bold"),
            axis.text.y.left = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15,face = "bold"),
            legend.position = "top")
    return(p)
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
        ggtitle(title1)+
        labs(color=gene)+
        ylab("Score")+
        scale_color_manual(values = c("red2", "blue4"))+
        geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    
    return(p)
  } 
  
  
}

immune_infiltration_CPTAC_rna_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_infiltration"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_infiltration"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(RNA)")+
    geom_boxplot(aes(x = immune_infiltration,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_CPTAC_protein_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_infiltration_protein"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_infiltration_protein"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(protein)")+
    geom_boxplot(aes(x = immune_infiltration_protein,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration_protein)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration_protein == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration_protein == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_CPTAC_rna_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_pathway"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_pathway"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (RNA)")+
    geom_boxplot(aes(x = immune_pathway,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_CPTAC_protein_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),wt_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,mut]
  wt_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,wt]
  
  mut_values = melt(mut_values)
  mut_values = mut_values[,-2]
  colnames(mut_values)[1] = "immune_pathway_protein"
  mut_values$groups = "Mutation"
  
  wt_values = melt(wt_values)
  wt_values = wt_values[,-2]
  colnames(wt_values)[1] = "immune_pathway_protein"
  wt_values$groups = "Wiltype"
  
  tmp_data = rbind(mut_values,wt_values)
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (protein)")+
    geom_boxplot(aes(x = immune_pathway_protein,y=value,fill=groups),alpha = 0.8)+
    theme_classic()+
    labs(fill = gene)+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway_protein)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway_protein == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway_protein == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

DEG_CPTAC_rna_pm = function(TCGA,gene,cancer_type,min.pct,Mut_type,Wild_type,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["rna"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["rna"]]),wt_patient_id)
  
  tmp_data = TCGA[[cancer_type]][["rna"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

DEG_CPTAC_protein_pm = function(TCGA,gene,cancer_type,min.pct,mut_patient_id,wt_patient_id){
  
  mut = intersect(colnames(TCGA[[cancer_type]][["protein"]]),mut_patient_id)
  wt = intersect(colnames(TCGA[[cancer_type]][["protein"]]),wt_patient_id)
  
  tmp_data = TCGA[[cancer_type]][["protein"]][,c(mut,wt)]
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("Mutation",length(mut)),rep("Wildtype",length(wt))),levels = c("Mutation","Wildtype"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("Mutation","Wildtype")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(Mutation-Wildtype,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
}

CPTAC_survival_pm = function(TCGA,cancer_type,gene,mut_patient_id,wt_patient_id){
  
  mut = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,mut_patient_id)
  wt = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,wt_patient_id)
  tmp_data = TCGA[[cancer_type]][["Survival"]]
  tmp_data$groups = NA
  tmp_data$groups[tmp_data$PatientID %in% mut] = "Mutation"
  tmp_data$groups[tmp_data$PatientID %in% wt] = "Wildtype"
  tmp_data = tmp_data[,-1]
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title=gene,
                  legend.labs = c("Mutation","Wildtype"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  p = (p1$plot/p1$cumevents)+plot_layout(design = layout)
  return(p)
  
}

#################################TCGA_Subtype########################################################
TCGA_cohort_cal_subtype = function(TCGA,Mutation_subtype,cancer_type,Mut_type,Wild_type){
  
  if("All" %in% Mut_type){
    mut_patient_id = unique(as.character(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode[TCGA[[cancer_type]][["maf"]]@data$Hugo_Symbol == Mutation_subtype]))
  }else{
    mut_patient_id = unique(as.character(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode[
      TCGA[[cancer_type]][["maf"]]@data$Hugo_Symbol == Mutation_subtype & 
        TCGA[[cancer_type]][["maf"]]@data$Variant_Classification %in% Mut_type
    ]))
  }
  
  if(Wild_type == "Others"){
    wt_patient_id = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode)),mut_patient_id) 
  }else{
    wt_patient_id = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode)),
                            unique(as.character(TCGA[[cancer_type]]$maf@data$Tumor_Sample_Barcode[TCGA[[cancer_type]]$maf@data$Hugo_Symbol == Mutation_subtype])))
  }
  
  return(list("mut" = mut_patient_id,"wt" = wt_patient_id))
}

immune_one_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,immune_module,selection,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
    pat_id = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),patient_id)
    pat_values = TCGA[[cancer_type]][[immune_module]][selection,pat_id]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
      cutoff = median(tmp_mean_exp)
      UP_name = pat_id[tmp_mean_exp > cutoff]
      DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
      cutoff = median(tmp_mean_exp)
      UP_name = pat_id[tmp_mean_exp > cutoff]
      DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    }
    
    tmp_data = data.frame(values=c(pat_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))))
    
  }else if(MT_OR_WT == "Mutation"){
    
    
    mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
      cutoff = median(tmp_mean_exp)
      UP_name = mut[tmp_mean_exp > cutoff]
      DOWN_name = mut[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
      cutoff = median(tmp_mean_exp)
      UP_name = mut[tmp_mean_exp > cutoff]
      DOWN_name = mut[tmp_mean_exp <= cutoff]
    }
    
    
    tmp_data = data.frame(values=c(mut_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))))
  }else{
    
    
    wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
      cutoff = median(tmp_mean_exp)
      UP_name = wt[tmp_mean_exp > cutoff]
      DOWN_name = wt[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
      cutoff = median(tmp_mean_exp)
      UP_name = wt[tmp_mean_exp > cutoff]
      DOWN_name = wt[tmp_mean_exp <= cutoff]
    }
    
    
    tmp_data = data.frame(values=c(wt_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))))
  }
  return(tmp_data)
}
immune_one_subtype = function(TCGA,immune_module,selection,outlier,tmp_data){
  
  res = wilcox.test(values~groups,data = tmp_data)
  if(grepl("protein",immune_module)){title1 = paste(selection,"(Protein)")}else{title1 = paste(selection,"(RNA)")}
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "High","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Low","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(title1)+
      labs(color="Genes expression")+
      ylab("Score")+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
            plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
            axis.text.x = element_text(size = 15,face = "bold"),
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 15,face = "bold"),
            axis.text.y.left = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15,face = "bold"),
            legend.position = "top")
    return(p)
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
        ggtitle(title1)+
        labs(color="Genes expression")+
        ylab("Score")+
        geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    
    return(p)
  } 
  
  
}

immune_one_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,immune_module,selection,mut_patient_id,wt_patient_id){
  
  
  patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
  pat_id = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),patient_id)
  pat_values = TCGA[[cancer_type]][[immune_module]][selection,pat_id]
  
  if(grepl("protein",immune_module)){
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  }else{
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  }
  
  All = data.frame(values=c(pat_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),module = "All")
  
  
  
  
  mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
  mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
  
  if(grepl("protein",immune_module)){
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
  }else{
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
  }
  
  
  Mutation = data.frame(values=c(mut_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),module = paste(Mutation_subtype,"Mutation"))
  
  
  
  wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
  wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
  
  if(grepl("protein",immune_module)){
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
  }else{
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
  }
  
  
  Wildtype = data.frame(values=c(wt_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),module = paste(Mutation_subtype,"Wildtype"))
  
  return(rbind(All,Mutation,Wildtype))
}
immune_one_subtype_compare = function(tmp_data_compare,selection,genes){
  tmp_data_compare$groups = paste(genes,tmp_data_compare$groups)
  p = ggplot(data = tmp_data_compare)+
    ggtitle(paste(selection,"(RNA)"))+
    geom_boxplot(aes(x = module,y = values,fill = groups),alpha = 0.8)+
    theme_classic()+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,size = 18,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
          axis.title.x= element_blank(),
          axis.title.y = element_text(size = 18,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.position = "right")
  n = 1
  for(i in unique(tmp_data_compare$module)){
    
    res = wilcox.test(values~groups,data = tmp_data_compare[tmp_data_compare$module == i,])
    max_val = max(tmp_data_compare[tmp_data_compare$module == i,"values"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 6)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_infiltration"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
  }
  return(tmp_data)
}
immune_infiltration_subtype = function(TCGA,tmp_data){
  
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(RNA)")+
    geom_boxplot(aes(x = immune_infiltration,y=value,fill=groups))+
    theme_classic()+
    labs(fill = "Genes expression")+
    ylab("Score")+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x= element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){
  
  patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
  pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),patient_id)
  pat_values = TCGA[[cancer_type]][["immune_infiltration"]][,pat_id]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
  cutoff = median(tmp_mean_exp)
  UP_name = pat_id[tmp_mean_exp > cutoff]
  DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  
  UP_values = melt(pat_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_infiltration"
  UP_values$groups = "High"
  
  DOWN_values = melt(pat_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_infiltration"
  DOWN_values$groups = "Low"
  
  All = rbind(UP_values,DOWN_values)
  All$module = "All"
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
  cutoff = median(tmp_mean_exp)
  UP_name = mut[tmp_mean_exp > cutoff]
  DOWN_name = mut[tmp_mean_exp <= cutoff]
  
  UP_values = melt(mut_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_infiltration"
  UP_values$groups = "High"
  
  DOWN_values = melt(mut_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_infiltration"
  DOWN_values$groups = "Low"
  
  Mutation = rbind(UP_values,DOWN_values)
  Mutation$module = paste(Mutation_subtype,"Mutation")
  
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
  wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
  cutoff = median(tmp_mean_exp)
  UP_name = wt[tmp_mean_exp > cutoff]
  DOWN_name = wt[tmp_mean_exp <= cutoff]
  
  UP_values = melt(wt_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_infiltration"
  UP_values$groups = "High"
  
  DOWN_values = melt(wt_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_infiltration"
  DOWN_values$groups = "Low"
  
  Wildtype = rbind(UP_values,DOWN_values)
  Wildtype$module = paste(Mutation_subtype,"Wildtype")
  
  return(rbind(All,Mutation,Wildtype))
}
immune_infiltration_subtype_compare = function(tmp_data_compare,genes){
  res = vector()
  for(i in unique(tmp_data_compare$module)){
    for(j in unique(tmp_data_compare$immune_infiltration)){
      tmp2 = tmp_data_compare[ tmp_data_compare$immune_infiltration == j & tmp_data_compare$module == i,]
      pvalue = wilcox.test(value~groups,data = tmp2)$p.value
      tmp3 = tmp2 %>% group_by(groups) %>% summarise(median(value))
      FC = as.numeric(tmp3[1,2] - tmp3[2,2]) 
      res = rbind(res,c(j,FC,pvalue,i))
    }
  }
  
  res = as.data.frame(res)
  colnames(res) = c("Immune","FC","pvalue","subtype")
  res$FC = as.numeric(res$FC)
  res$pvalue = as.numeric(res$pvalue)
  res$Significance = ifelse(res$pvalue <= 0.05,"P.value <= 0.05","P.value > 0.05")
  
  p = ggplot(data = res)+
    ggtitle("Infiltrating immune cells(RNA)")+
    theme_bw()+
    geom_point(aes(x = Immune,y = subtype,size = -log10(pvalue),fill = FC,color = Significance),shape = 21,)+
    scale_color_manual(values = c("P.value <= 0.05" = "red","P.value > 0.05" = "grey"))+
    # scale_fill_viridis(limits = c(-max(res$FC),max(res$FC)),breaks = c(-max(res$FC),max(res$FC)),labels = c(paste(genes,"Low ↑↑"),paste(genes,"High ↑↑")))+
    scale_fill_gradient2(low = "blue3",mid = "white", high = "yellow3",limits = c(-max(abs(res$FC)),max(abs(res$FC))),breaks = c(-max(abs(res$FC)),max(abs(res$FC))),labels = c(paste(genes,"Low (Up)"),paste(genes,"High (Up)")))+
    theme(axis.text.x = element_text(angle = 30,size = 12,hjust = 1,face = "bold"),plot.title = element_text(size = 22,face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.title = element_blank(),
          axis.text.y = element_text(face = "bold",size = 15))
  
  return(p)
  
  
}

immune_signature_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_pathway"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
  }
}
immune_signature_subtype = function(TCGA,tmp_data){
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (RNA)")+
    geom_boxplot(aes(x = immune_pathway,y=value,fill=groups))+
    theme_classic()+
    labs(fill = "Gene expression")+
    ylab("Score")+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){
  
  patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
  pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),patient_id)
  pat_values = TCGA[[cancer_type]][["immune_pathway"]][,pat_id]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
  cutoff = median(tmp_mean_exp)
  UP_name = pat_id[tmp_mean_exp > cutoff]
  DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  
  UP_values = melt(pat_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway"
  UP_values$groups = "High"
  
  DOWN_values = melt(pat_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway"
  DOWN_values$groups = "Low"
  
  All = rbind(UP_values,DOWN_values)
  All$module = "All"
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
  cutoff = median(tmp_mean_exp)
  UP_name = mut[tmp_mean_exp > cutoff]
  DOWN_name = mut[tmp_mean_exp <= cutoff]
  
  UP_values = melt(mut_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway"
  UP_values$groups = "High"
  
  DOWN_values = melt(mut_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway"
  DOWN_values$groups = "Low"
  
  Mutation = rbind(UP_values,DOWN_values)
  Mutation$module = paste(Mutation_subtype,"Mutation")
  
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
  wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
  cutoff = median(tmp_mean_exp)
  UP_name = wt[tmp_mean_exp > cutoff]
  DOWN_name = wt[tmp_mean_exp <= cutoff]
  
  UP_values = melt(wt_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway"
  UP_values$groups = "High"
  
  DOWN_values = melt(wt_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway"
  DOWN_values$groups = "Low"
  
  Wildtype = rbind(UP_values,DOWN_values)
  Wildtype$module = paste(Mutation_subtype,"Wildtype")
  
  return(rbind(All,Mutation,Wildtype))
}
immune_signature_subtype_compare = function(tmp_data_compare,genes){
  res = vector()
  for(i in unique(tmp_data_compare$module)){
    for(j in unique(tmp_data_compare$immune_pathway)){
      tmp2 = tmp_data_compare[ tmp_data_compare$immune_pathway == j & tmp_data_compare$module == i,]
      pvalue = wilcox.test(value~groups,data = tmp2)$p.value
      tmp3 = tmp2 %>% group_by(groups) %>% summarise(median(value))
      FC = as.numeric(tmp3[1,2] - tmp3[2,2]) 
      res = rbind(res,c(j,FC,pvalue,i))
    }
  }
  
  res = as.data.frame(res)
  colnames(res) = c("Immune","FC","pvalue","subtype")
  res$FC = as.numeric(res$FC)
  res$pvalue = as.numeric(res$pvalue)
  res$Significance = ifelse(res$pvalue <= 0.05,"P.value <= 0.05","P.value > 0.05")
  
  p = ggplot(data = res)+
    ggtitle("Infiltrating immune cells(RNA)")+
    theme_bw()+
    geom_point(aes(x = Immune,y = subtype,size = -log10(pvalue),fill = FC,color = Significance),shape = 21,)+
    scale_color_manual(values = c("P.value <= 0.05" = "red","P.value > 0.05" = "grey"))+
    # scale_fill_viridis(limits = c(-max(res$FC),max(res$FC)),breaks = c(-max(res$FC),max(res$FC)),labels = c(paste(genes,"Low ↑↑"),paste(genes,"High ↑↑")))+
    scale_fill_gradient2(low = "blue3",mid = "white", high = "yellow3",limits = c(-max(abs(res$FC)),max(abs(res$FC))),breaks = c(-max(abs(res$FC)),max(abs(res$FC))),labels = c(paste(genes,"Low (Up)"),paste(genes,"High (Up)")))+
    theme(axis.text.x = element_text(angle = 30,size = 12,hjust = 1,face = "bold"),plot.title = element_text(size = 22,face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.title = element_blank(),
          axis.text.y = element_text(face = "bold",size = 15))
  
  return(p)
  
  
}

DEG_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["rna"]]),patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    
    tmp_data = TCGA[[cancer_type]][["rna"]][,c(UP_name,DOWN_name)]
    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][["rna"]]),mut_patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    
    tmp_data = TCGA[[cancer_type]][["rna"]][,c(UP_name,DOWN_name)]
    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["rna"]]),wt_patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    tmp_data = TCGA[[cancer_type]][["rna"]][,c(UP_name,DOWN_name)]
    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
  }
}
DEG_subtype = function(TCGA,min.pct,tmp_data){
  UP_name = tmp_data$UP
  DOWN_name = tmp_data$DOWN
  tmp_data = tmp_data$tmp_data
  
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),levels = c("High","Low"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("High","Low")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(High-Low,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
  
}

TCGA_survival_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
    pat_id = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),patient_id)
    pat_id2 = intersect(pat_id,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["clinical_detial"]][pat_id2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id2])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id2[tmp_mean_exp > cutoff]
    DOWN_name = pat_id2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[UP_name,"groups"] = "High"
    tmp_data[DOWN_name,"groups"] = "low"
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),mut_patient_id)
    mut2 = intersect(mut,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["clinical_detial"]][mut2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut2])
    cutoff = median(tmp_mean_exp)
    UP_name = mut2[tmp_mean_exp > cutoff]
    DOWN_name = mut2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[UP_name,"groups"] = "High"
    tmp_data[DOWN_name,"groups"] = "low"
    
  }else{
    
    wt = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),wt_patient_id)
    wt2 = intersect(wt,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["clinical_detial"]][wt2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt2])
    cutoff = median(tmp_mean_exp)
    UP_name = wt2[tmp_mean_exp > cutoff]
    DOWN_name = wt2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[UP_name,"groups"] = "High"
    tmp_data[DOWN_name,"groups"] = "low"
  }
  
  tmp_data = tmp_data[,-1]
  return(tmp_data)
}
TCGA_survival_subtype = function(TCGA,cancer_type,tmp_data){
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title="Genes expression",
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  
  
  
  if(cancer_type != "LAML"){
    
    fit <- do.call(survfit, list(Surv(DSS.time,DSS)~tmp_data[,"groups"],data = tmp_data))
    p2 = ggsurvplot(fit,data = tmp_data,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Disease-specific survival (DSS)",
                    legend.title="Genes expression",
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    
    fit <- do.call(survfit, list(Surv(PFI.time,PFI)~tmp_data[,"groups"],data = tmp_data))
    p3 = ggsurvplot(fit,data = tmp_data,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Progression-free interval (PFI)",
                    legend.title="Genes expression",
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  }

  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  if(cancer_type != "LAML"){
    p = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)  | (p3$plot/p3$cumevents)+plot_layout(design = layout)
    return(p)
  }else{
    p = (p1$plot/p1$cumevents)+plot_layout(design = layout) 
    return(p)
  }  
  
}

TCGA_survival_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){

    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$Tumor_Sample_Barcode)
    pat_id = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),patient_id)
    pat_id2 = intersect(pat_id,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["clinical_detial"]][pat_id2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id2])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id2[tmp_mean_exp > cutoff]
    DOWN_name = pat_id2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[UP_name,"groups"] = "High"
    tmp_data[DOWN_name,"groups"] = "low"
    All = tmp_data
    All$module = "All"

    
    mut = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),mut_patient_id)
    mut2 = intersect(mut,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["clinical_detial"]][mut2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut2])
    cutoff = median(tmp_mean_exp)
    UP_name = mut2[tmp_mean_exp > cutoff]
    DOWN_name = mut2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[UP_name,"groups"] = "High"
    tmp_data[DOWN_name,"groups"] = "low"
    Mutation = tmp_data
    Mutation$module = paste(Mutation_subtype,"Mutation")

    
    wt = intersect(rownames(TCGA[[cancer_type]][["clinical_detial"]]),wt_patient_id)
    wt2 = intersect(wt,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["clinical_detial"]][wt2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt2])
    cutoff = median(tmp_mean_exp)
    UP_name = wt2[tmp_mean_exp > cutoff]
    DOWN_name = wt2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[UP_name,"groups"] = "High"
    tmp_data[DOWN_name,"groups"] = "low"
    Wildtype = tmp_data
    Wildtype$module = paste(Mutation_subtype,"Wildtype")
  
    tmp_data = rbind(All,Mutation,Wildtype)
    tmp_data = tmp_data[,-1]
    return(tmp_data)
}
TCGA_survival_subtype_compare = function(tmp_data_compare,cancer_type,genes){
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  
  tmp2 = tmp_data_compare[ tmp_data_compare$module == "All",]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p1 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)\n(All)",
                  legend.title=genes,
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  tmp2 = tmp_data_compare[ grepl(pattern = "Mutation",tmp_data_compare$module),]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p2 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title = paste("Overall Survival (OS)","\n(",tmp2$module[1],")",sep = ""),
                  legend.title=genes,
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  tmp2 = tmp_data_compare[ grepl(pattern = "Wildtype",tmp_data_compare$module),]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p3 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title = paste("Overall Survival (OS)","\n(",tmp2$module[1],")",sep = ""),
                  legend.title=genes,
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  pOS = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)  | (p3$plot/p3$cumevents)+plot_layout(design = layout)
  
  if(cancer_type != "LAML"){
    tmp2 = tmp_data_compare[ grepl(pattern = "All",tmp_data_compare$module),]
    fit <- do.call(survfit, list(Surv(DSS.time,DSS)~tmp2[,"groups"],data = tmp2))
    p1 = ggsurvplot(fit,data = tmp2,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Disease-specific survival (DSS)\n(All)",
                    legend.title=genes,
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    tmp2 = tmp_data_compare[ grepl(pattern = "Mutation",tmp_data_compare$module),]
    fit <- do.call(survfit, list(Surv(DSS.time,DSS)~tmp2[,"groups"],data = tmp2))
    p2 = ggsurvplot(fit,data = tmp2,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title=paste("Disease-specific survival (DSS)","\n(",tmp2$module[1],")",sep = ""),
                    legend.title=genes,
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    tmp2 = tmp_data_compare[ grepl(pattern = "Wildtype",tmp_data_compare$module),]
    fit <- do.call(survfit, list(Surv(DSS.time,DSS)~tmp2[,"groups"],data = tmp2))
    p3 = ggsurvplot(fit,data = tmp2,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title=paste("Disease-specific survival (DSS)","\n(",tmp2$module[1],")",sep = ""),
                    legend.title=genes,
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    pDSS = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)  | (p3$plot/p3$cumevents)+plot_layout(design = layout)
    
    
    tmp2 = tmp_data_compare[ grepl(pattern = "All",tmp_data_compare$module),]
    fit <- do.call(survfit, list(Surv(PFI.time,PFI)~tmp2[,"groups"],data = tmp2))
    p1 = ggsurvplot(fit,data = tmp2,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title="Progression-free interval (PFI)\n(All)",
                    legend.title=genes,
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    tmp2 = tmp_data_compare[ grepl(pattern = "Mutation",tmp_data_compare$module),]
    fit <- do.call(survfit, list(Surv(PFI.time,PFI)~tmp2[,"groups"],data = tmp2))
    p2 = ggsurvplot(fit,data = tmp2,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title=paste("Progression-free interval (PFI)","\n(",tmp2$module[1],")",sep = ""),
                    legend.title=genes,
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    tmp2 = tmp_data_compare[ grepl(pattern = "Wildtype",tmp_data_compare$module),]
    fit <- do.call(survfit, list(Surv(PFI.time,PFI)~tmp2[,"groups"],data = tmp2))
    p3 = ggsurvplot(fit,data = tmp2,
                    xlab = "day",
                    pval = T,
                    risk.table = F,
                    cumevents = T,
                    palette = c("red2", "blue4"),
                    title=paste("Progression-free interval (PFI)","\n(",tmp2$module[1],")",sep = ""),
                    legend.title=genes,
                    legend.labs = c("High","low"),
                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
    
    pPFI = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout)  | (p3$plot/p3$cumevents)+plot_layout(design = layout)
    
  }
  
  if(cancer_type != "LAML"){
    p = pOS / pDSS /pPFI
    return(p)
  }else{
    return(pOS)
  } 
  
}
#################################CPTAC_Subtype########################################################
CPTAC_cohort_cal_subtype = function(TCGA,Mutation_subtype,cancer_type,Mut_type,Wild_type){
  
  if("All" %in% Mut_type){
    mut_patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA[TCGA[[cancer_type]][["maf"]]@data$Hugo_Symbol == Mutation_subtype])
  }else{
    mut_patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA[
      TCGA[[cancer_type]][["maf"]]@data$Hugo_Symbol == Mutation_subtype & 
        TCGA[[cancer_type]][["maf"]]@data$Variant_Classification %in% Mut_type
    ])
  }
  
  if(Wild_type == "Others"){
    wt_patient_id = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA)),mut_patient_id) 
  }else{
    wt_patient_id = setdiff(x = unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA)),
                            unique(as.character(TCGA[[cancer_type]]$maf@data$TCGA[TCGA[[cancer_type]]$maf@data$Hugo_Symbol == Mutation_subtype])))
  }
  
  return(list("mut" = mut_patient_id,"wt" = wt_patient_id))
}

immune_one_CPTAC_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,immune_module,selection,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),patient_id)
    pat_values = TCGA[[cancer_type]][[immune_module]][selection,pat_id]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
      cutoff = median(tmp_mean_exp)
      UP_name = pat_id[tmp_mean_exp > cutoff]
      DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
      cutoff = median(tmp_mean_exp)
      UP_name = pat_id[tmp_mean_exp > cutoff]
      DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    }
    
    
    tmp_data = data.frame(values=c(pat_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name)))) 
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
      cutoff = median(tmp_mean_exp)
      UP_name = mut[tmp_mean_exp > cutoff]
      DOWN_name = mut[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
      cutoff = median(tmp_mean_exp)
      UP_name = mut[tmp_mean_exp > cutoff]
      DOWN_name = mut[tmp_mean_exp <= cutoff]
    }
    
    
    tmp_data = data.frame(values=c(mut_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))))
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
      cutoff = median(tmp_mean_exp)
      UP_name = wt[tmp_mean_exp > cutoff]
      DOWN_name = wt[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
      cutoff = median(tmp_mean_exp)
      UP_name = wt[tmp_mean_exp > cutoff]
      DOWN_name = wt[tmp_mean_exp <= cutoff]
    }
    
    
    tmp_data = data.frame(values=c(wt_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))))
  }
  
  return(tmp_data)
  
}
immune_one_CPTAC_subtype = function(TCGA,immune_module,selection,outlier,tmp_data){
  
 
  res = wilcox.test(values~groups,data = tmp_data)
  if(grepl("protein",immune_module)){title1 = paste(selection,"(Protein)")}else{title1 = paste(selection,"(RNA)")}
  if(outlier){
    
    lim1 = boxplot.stats(tmp_data[tmp_data$groups == "High","values"])$stats[c(1, 5)]
    lim2 = boxplot.stats(tmp_data[tmp_data$groups == "Low","values"])$stats[c(1, 5)]
    ylim_all = c(lim1,lim2)
    ylim1 = c(min(ylim_all),max(ylim_all))
    
    p0 = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.color = NA)+
      ggtitle(title1)+
      labs(color="Genes expression")+
      ylab("Score")+
      geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max(ylim1)-min(ylim1))/100))+
      theme_classic()
    
    
    p1 = p0 + 
      coord_cartesian(ylim = c(ylim1[1],ylim1[2]*1.1))+
      annotate("text",x = 1.5,y = max(ylim1*1.1),label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)
    
    
    p = p1+ theme(
            plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
            axis.text.x = element_text(size = 15,face = "bold"),
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 15,face = "bold"),
            axis.text.y.left = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15,face = "bold"),
            legend.position = "top")
    return(p)
  }else{
    
    max_val = max(tmp_data$values)
    min_val = min(tmp_data$values)
    p = ggplot(data = tmp_data)+geom_boxplot(mapping = aes(x = groups,y = values,color = groups),alpha = 0.5,outlier.colour = NA)+
        ggtitle(title1)+
        labs(color="Genes expression")+
        ylab("Score")+
        geom_point(mapping = aes(x = groups,y = values,color = groups),position = position_jitter(width = 0.2,height = (max_val-min_val)/100))+
        theme_classic()+
        annotate("text",x = 1.5,y = max_val,label = paste("Wilcox : P value = ",round(res$p.value,digits = 4)),size = 6)+
        theme(
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.text.x = element_text(size = 15,face = "bold"),
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y.left = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
    return(p)
    
  } 
  
  
}

immune_one_CPTAC_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,immune_module,selection,mut_patient_id,wt_patient_id){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),patient_id)
    pat_values = TCGA[[cancer_type]][[immune_module]][selection,pat_id]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
      cutoff = median(tmp_mean_exp)
      UP_name = pat_id[tmp_mean_exp > cutoff]
      DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
      cutoff = median(tmp_mean_exp)
      UP_name = pat_id[tmp_mean_exp > cutoff]
      DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    }
    
    
    All = data.frame(values=c(pat_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),module = "All")
    
    mut = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][[immune_module]][selection,mut]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
      cutoff = median(tmp_mean_exp)
      UP_name = mut[tmp_mean_exp > cutoff]
      DOWN_name = mut[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
      cutoff = median(tmp_mean_exp)
      UP_name = mut[tmp_mean_exp > cutoff]
      DOWN_name = mut[tmp_mean_exp <= cutoff]
    }
    
    
    Mutation = data.frame(values=c(mut_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),module = paste(Mutation_subtype,"Mutation"))
    
    wt = intersect(colnames(TCGA[[cancer_type]][[immune_module]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][[immune_module]][selection,wt]
    
    if(grepl("protein",immune_module)){
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
      cutoff = median(tmp_mean_exp)
      UP_name = wt[tmp_mean_exp > cutoff]
      DOWN_name = wt[tmp_mean_exp <= cutoff]
    }else{
      tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
      cutoff = median(tmp_mean_exp)
      UP_name = wt[tmp_mean_exp > cutoff]
      DOWN_name = wt[tmp_mean_exp <= cutoff]
    }
    
    
    Wildtype = data.frame(values=c(wt_values[c(UP_name,DOWN_name)]),groups = c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),module = paste(Mutation_subtype,"Wildtype"))

  
  return(rbind(All,Mutation,Wildtype))
  
}
immune_one_CPTAC_subtype_compare = function(tmp_data_compare,immune_module,selection,genes){
  
  if(grepl("protein",immune_module)){title1 = paste(selection,"(Protein)")}else{title1 = paste(selection,"(RNA)")}
  tmp_data_compare$groups = paste(genes,tmp_data_compare$groups)
  p = ggplot(data = tmp_data_compare)+
    ggtitle(title1)+
    geom_boxplot(aes(x = module,y = values,fill = groups),alpha = 0.8)+
    theme_classic()+
    ylab("Score")+
    scale_fill_manual(values = c("red2", "blue4"))+
    theme(axis.text.x = element_text(hjust = 1,size = 18,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 22,face = "bold"),
          axis.title.x= element_blank(),
          axis.title.y = element_text(size = 18,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.position = "right")
  n = 1
  for(i in unique(tmp_data_compare$module)){
    
    res = wilcox.test(values~groups,data = tmp_data_compare[tmp_data_compare$module == i,])
    max_val = max(tmp_data_compare[tmp_data_compare$module == i,"values"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 6)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_CPTAC_rna_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_infiltration"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
  }
  return(tmp_data)
}
immune_infiltration_CPTAC_rna_subtype = function(TCGA,tmp_data){
  
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(RNA)")+
    geom_boxplot(aes(x = immune_infiltration,y=value,fill=groups))+
    theme_classic()+
    labs(fill = "Genes expression")+
    ylab("Score")+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_CPTAC_rna_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_infiltration"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    All = rbind(UP_values,DOWN_values)
    All$module = "All"

    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_infiltration"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    Mutation = rbind(UP_values,DOWN_values)
    Mutation$module = paste(Mutation_subtype,"Mutation")

    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_infiltration"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration"
    DOWN_values$groups = "Low"
    
    Wildtype = rbind(UP_values,DOWN_values)
    Wildtype$module = paste(Mutation_subtype,"Wildtype")
    
  return(rbind(All,Mutation,Wildtype))
}
immune_infiltration_CPTAC_rna_subtype_compare = function(tmp_data_compare,genes){
  res = vector()
  for(i in unique(tmp_data_compare$module)){
    for(j in unique(tmp_data_compare$immune_infiltration)){
      tmp2 = tmp_data_compare[ tmp_data_compare$immune_infiltration == j & tmp_data_compare$module == i,]
      pvalue = wilcox.test(value~groups,data = tmp2)$p.value
      tmp3 = tmp2 %>% group_by(groups) %>% summarise(median(value))
      FC = as.numeric(tmp3[1,2] - tmp3[2,2]) 
      res = rbind(res,c(j,FC,pvalue,i))
    }
  }
  
  res = as.data.frame(res)
  colnames(res) = c("Immune","FC","pvalue","subtype")
  res$FC = as.numeric(res$FC)
  res$pvalue = as.numeric(res$pvalue)
  res$Significance = ifelse(res$pvalue <= 0.05,"P.value <= 0.05","P.value > 0.05")
  
  p = ggplot(data = res)+
    ggtitle("Infiltrating immune cells(RNA)")+
    theme_bw()+
    geom_point(aes(x = Immune,y = subtype,size = -log10(pvalue),fill = FC,color = Significance),shape = 21,)+
    scale_color_manual(values = c("P.value <= 0.05" = "red","P.value > 0.05" = "grey"))+
    # scale_fill_viridis(limits = c(-max(res$FC),max(res$FC)),breaks = c(-max(res$FC),max(res$FC)),labels = c(paste(genes,"Low ↑↑"),paste(genes,"High ↑↑")))+
    scale_fill_gradient2(low = "blue3",mid = "white", high = "yellow3",limits = c(-max(abs(res$FC)),max(abs(res$FC))),breaks = c(-max(abs(res$FC)),max(abs(res$FC))),labels = c(paste(genes,"Low (Up)"),paste(genes,"High (Up)")))+
    theme(axis.text.x = element_text(angle = 30,size = 12,hjust = 1,face = "bold"),plot.title = element_text(size = 22,face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.title = element_blank(),
          axis.text.y = element_text(face = "bold",size = 15))
  
  return(p)
  
  
}

immune_infiltration_CPTAC_protein_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration_protein"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration_protein"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else if(MT_OR_WT == "Mutation"){
    
    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration_protein"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration_protein"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_infiltration_protein"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_infiltration_protein"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
  }
  
  return(tmp_data)
}
immune_infiltration_CPTAC_protein_subtype = function(TCGA,tmp_data){
 
  p = ggplot(data = tmp_data)+
    ggtitle("Infiltrating immune cells(Protein)")+
    geom_boxplot(aes(x = immune_infiltration_protein,y=value,fill=groups))+
    theme_classic()+
    labs(fill = "Genes expression")+
    ylab("Score")+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 12,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_infiltration_protein)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_infiltration_protein == i,])
    max_val = max(tmp_data[tmp_data$immune_infiltration_protein == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_infiltration_CPTAC_protein_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){
  
  patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
  pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),patient_id)
  pat_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,pat_id]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
  cutoff = median(tmp_mean_exp)
  UP_name = pat_id[tmp_mean_exp > cutoff]
  DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  
  UP_values = melt(pat_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_infiltration_protein"
  UP_values$groups = "High"
  
  DOWN_values = melt(pat_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_infiltration_protein"
  DOWN_values$groups = "Low"
  
  All = rbind(UP_values,DOWN_values)
  All$module = "All"
  
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),mut_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,mut]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
  cutoff = median(tmp_mean_exp)
  UP_name = mut[tmp_mean_exp > cutoff]
  DOWN_name = mut[tmp_mean_exp <= cutoff]
  
  UP_values = melt(mut_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_infiltration_protein"
  UP_values$groups = "High"
  
  DOWN_values = melt(mut_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_infiltration_protein"
  DOWN_values$groups = "Low"
  
  Mutation = rbind(UP_values,DOWN_values)
  Mutation$module = paste(Mutation_subtype,"Mutation")
  
  
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_infiltration_protein"]]),wt_patient_id)
  wt_values = TCGA[[cancer_type]][["immune_infiltration_protein"]][,wt]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
  cutoff = median(tmp_mean_exp)
  UP_name = wt[tmp_mean_exp > cutoff]
  DOWN_name = wt[tmp_mean_exp <= cutoff]
  
  UP_values = melt(wt_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_infiltration_protein"
  UP_values$groups = "High"
  
  DOWN_values = melt(wt_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_infiltration_protein"
  DOWN_values$groups = "Low"
  
  Wildtype = rbind(UP_values,DOWN_values)
  Wildtype$module = paste(Mutation_subtype,"Wildtype")
  
  return(rbind(All,Mutation,Wildtype))
}
immune_infiltration_CPTAC_protein_subtype_compare = function(tmp_data_compare,genes){
  res = vector()
  for(i in unique(tmp_data_compare$module)){
    for(j in unique(tmp_data_compare$immune_infiltration_protein)){
      tmp2 = tmp_data_compare[ tmp_data_compare$immune_infiltration_protein == j & tmp_data_compare$module == i,]
      pvalue = wilcox.test(value~groups,data = tmp2)$p.value
      tmp3 = tmp2 %>% group_by(groups) %>% summarise(median(value))
      FC = as.numeric(tmp3[1,2] - tmp3[2,2]) 
      res = rbind(res,c(j,FC,pvalue,i))
    }
  }
  
  res = as.data.frame(res)
  colnames(res) = c("Immune","FC","pvalue","subtype")
  res$FC = as.numeric(res$FC)
  res$pvalue = as.numeric(res$pvalue)
  res$Significance = ifelse(res$pvalue <= 0.05,"P.value <= 0.05","P.value > 0.05")
  
  p = ggplot(data = res)+
    ggtitle("Infiltrating immune cells(Protein)")+
    theme_bw()+
    geom_point(aes(x = Immune,y = subtype,size = -log10(pvalue),fill = FC,color = Significance),shape = 21,)+
    scale_color_manual(values = c("P.value <= 0.05" = "red","P.value > 0.05" = "grey"))+
    # scale_fill_viridis(limits = c(-max(res$FC),max(res$FC)),breaks = c(-max(res$FC),max(res$FC)),labels = c(paste(genes,"Low ↑↑"),paste(genes,"High ↑↑")))+
    scale_fill_gradient2(low = "blue3",mid = "white", high = "yellow3",limits = c(-max(abs(res$FC)),max(abs(res$FC))),breaks = c(-max(abs(res$FC)),max(abs(res$FC))),labels = c(paste(genes,"Low (Up)"),paste(genes,"High (Up)")))+
    theme(axis.text.x = element_text(angle = 30,size = 12,hjust = 1,face = "bold"),plot.title = element_text(size = 22,face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.title = element_blank(),
          axis.text.y = element_text(face = "bold",size = 15))
  
  return(p)
  
  
}

immune_signature_CPTAC_rna_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_pathway"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else if(MT_OR_WT == "Mutation"){
    
    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
  }
  return(tmp_data)
}
immune_signature_CPTAC_rna_subtype = function(TCGA,tmp_data){
  
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (RNA)")+
    geom_boxplot(aes(x = immune_pathway,y=value,fill=groups))+
    theme_classic()+
    labs(fill = "Gene expression")+
    ylab("Score")+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_CPTAC_rna_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){
  
  patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
  pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),patient_id)
  pat_values = TCGA[[cancer_type]][["immune_pathway"]][,pat_id]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
  cutoff = median(tmp_mean_exp)
  UP_name = pat_id[tmp_mean_exp > cutoff]
  DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  
  UP_values = melt(pat_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway"
  UP_values$groups = "High"
  
  DOWN_values = melt(pat_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway"
  DOWN_values$groups = "Low"
  
  All = rbind(UP_values,DOWN_values)
  All$module = "All"
  
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),mut_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway"]][,mut]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
  cutoff = median(tmp_mean_exp)
  UP_name = mut[tmp_mean_exp > cutoff]
  DOWN_name = mut[tmp_mean_exp <= cutoff]
  
  UP_values = melt(mut_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway"
  UP_values$groups = "High"
  
  DOWN_values = melt(mut_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway"
  DOWN_values$groups = "Low"
  
  Mutation = rbind(UP_values,DOWN_values)
  Mutation$module = paste(Mutation_subtype,"Mutation")
  
  
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway"]]),wt_patient_id)
  wt_values = TCGA[[cancer_type]][["immune_pathway"]][,wt]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
  cutoff = median(tmp_mean_exp)
  UP_name = wt[tmp_mean_exp > cutoff]
  DOWN_name = wt[tmp_mean_exp <= cutoff]
  
  UP_values = melt(wt_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway"
  UP_values$groups = "High"
  
  DOWN_values = melt(wt_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway"
  DOWN_values$groups = "Low"
  
  Wildtype = rbind(UP_values,DOWN_values)
  Wildtype$module = paste(Mutation_subtype,"Wildtype")
  
  return(rbind(All,Mutation,Wildtype))
}
immune_signature_CPTAC_rna_subtype_compare = function(tmp_data_compare,genes){
  res = vector()
  for(i in unique(tmp_data_compare$module)){
    for(j in unique(tmp_data_compare$immune_pathway)){
      tmp2 = tmp_data_compare[ tmp_data_compare$immune_pathway == j & tmp_data_compare$module == i,]
      pvalue = wilcox.test(value~groups,data = tmp2)$p.value
      tmp3 = tmp2 %>% group_by(groups) %>% summarise(median(value))
      FC = as.numeric(tmp3[1,2] - tmp3[2,2]) 
      res = rbind(res,c(j,FC,pvalue,i))
    }
  }
  
  res = as.data.frame(res)
  colnames(res) = c("Immune","FC","pvalue","subtype")
  res$FC = as.numeric(res$FC)
  res$pvalue = as.numeric(res$pvalue)
  res$Significance = ifelse(res$pvalue <= 0.05,"P.value <= 0.05","P.value > 0.05")
  
  p = ggplot(data = res)+
    ggtitle("Immune-related signatures (RNA)")+
    theme_bw()+
    geom_point(aes(x = Immune,y = subtype,size = -log10(pvalue),fill = FC,color = Significance),shape = 21,)+
    scale_color_manual(values = c("P.value <= 0.05" = "red","P.value > 0.05" = "grey"))+
    # scale_fill_viridis(limits = c(-max(res$FC),max(res$FC)),breaks = c(-max(res$FC),max(res$FC)),labels = c(paste(genes,"Low ↑↑"),paste(genes,"High ↑↑")))+
    scale_fill_gradient2(low = "blue3",mid = "white", high = "yellow3",limits = c(-max(abs(res$FC)),max(abs(res$FC))),breaks = c(-max(abs(res$FC)),max(abs(res$FC))),labels = c(paste(genes,"Low (Up)"),paste(genes,"High (Up)")))+
    theme(axis.text.x = element_text(angle = 30,size = 12,hjust = 1,face = "bold"),plot.title = element_text(size = 22,face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.title = element_blank(),
          axis.text.y = element_text(face = "bold",size = 15))
  
  return(p)
  
  
}

immune_signature_CPTAC_protein_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),patient_id)
    pat_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,pat_id]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    UP_values = melt(pat_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway_protein"
    UP_values$groups = "High"
    
    DOWN_values = melt(pat_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway_protein"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else if(MT_OR_WT == "Mutation"){
    
    
    mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),mut_patient_id)
    mut_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,mut]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    UP_values = melt(mut_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway_protein"
    UP_values$groups = "High"
    
    DOWN_values = melt(mut_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway_protein"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
    
  }else{
    
    
    wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),wt_patient_id)
    wt_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,wt]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    UP_values = melt(wt_values[,UP_name])
    UP_values = UP_values[,-2]
    colnames(UP_values)[1] = "immune_pathway_protein"
    UP_values$groups = "High"
    
    DOWN_values = melt(wt_values[,DOWN_name])
    DOWN_values = DOWN_values[,-2]
    colnames(DOWN_values)[1] = "immune_pathway_protein"
    DOWN_values$groups = "Low"
    
    tmp_data = rbind(UP_values,DOWN_values)
  }
  return(tmp_data)
}
immune_signature_CPTAC_protein_subtype = function(TCGA,tmp_data){
 
  p = ggplot(data = tmp_data)+
    ggtitle("Immune-related signatures (Protein)")+
    geom_boxplot(aes(x = immune_pathway_protein,y=value,fill=groups))+
    theme_classic()+
    labs(fill = "Gene expression")+
    ylab("Score")+
    theme(axis.text.x = element_text(hjust = 1,angle = 15,size = 15,face = "bold"),
          plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15,face = "bold"),
          legend.position = "top")
  n = 1
  for(i in unique(tmp_data$immune_pathway_protein)){
    
    res = wilcox.test(value~groups,data = tmp_data[tmp_data$immune_pathway_protein == i,])
    max_val = max(tmp_data[tmp_data$immune_pathway_protein == i,"value"])
    p = p+
      annotate("text",x = n,y = max_val*1.05,label = paste("Pvalue=",round(res$p.value,digits = 4),sep = ""),size = 3)
    n = n + 1
  }
  
  return(p)
}

immune_signature_CPTAC_protein_subtype_tmpdata_compare = function(TCGA,cancer_type,Mutation_subtype,genes,mut_patient_id,wt_patient_id){
  
  patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
  pat_id = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),patient_id)
  pat_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,pat_id]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
  cutoff = median(tmp_mean_exp)
  UP_name = pat_id[tmp_mean_exp > cutoff]
  DOWN_name = pat_id[tmp_mean_exp <= cutoff]
  
  UP_values = melt(pat_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway_protein"
  UP_values$groups = "High"
  
  DOWN_values = melt(pat_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway_protein"
  DOWN_values$groups = "Low"
  
  All = rbind(UP_values,DOWN_values)
  All$module = "All"
  
  
  mut = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),mut_patient_id)
  mut_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,mut]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
  cutoff = median(tmp_mean_exp)
  UP_name = mut[tmp_mean_exp > cutoff]
  DOWN_name = mut[tmp_mean_exp <= cutoff]
  
  UP_values = melt(mut_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway_protein"
  UP_values$groups = "High"
  
  DOWN_values = melt(mut_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway_protein"
  DOWN_values$groups = "Low"
  
  Mutation = rbind(UP_values,DOWN_values)
  Mutation$module = paste(Mutation_subtype,"Mutation")
  
  
  wt = intersect(colnames(TCGA[[cancer_type]][["immune_pathway_protein"]]),wt_patient_id)
  wt_values = TCGA[[cancer_type]][["immune_pathway_protein"]][,wt]
  
  tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
  cutoff = median(tmp_mean_exp)
  UP_name = wt[tmp_mean_exp > cutoff]
  DOWN_name = wt[tmp_mean_exp <= cutoff]
  
  UP_values = melt(wt_values[,UP_name])
  UP_values = UP_values[,-2]
  colnames(UP_values)[1] = "immune_pathway_protein"
  UP_values$groups = "High"
  
  DOWN_values = melt(wt_values[,DOWN_name])
  DOWN_values = DOWN_values[,-2]
  colnames(DOWN_values)[1] = "immune_pathway_protein"
  DOWN_values$groups = "Low"
  
  Wildtype = rbind(UP_values,DOWN_values)
  Wildtype$module = paste(Mutation_subtype,"Wildtype")
  
  return(rbind(All,Mutation,Wildtype))
}
immune_signature_CPTAC_protein_subtype_compare = function(tmp_data_compare,genes){
  res = vector()
  for(i in unique(tmp_data_compare$module)){
    for(j in unique(tmp_data_compare$immune_pathway_protein)){
      tmp2 = tmp_data_compare[ tmp_data_compare$immune_pathway_protein == j & tmp_data_compare$module == i,]
      pvalue = wilcox.test(value~groups,data = tmp2)$p.value
      tmp3 = tmp2 %>% group_by(groups) %>% summarise(median(value))
      FC = as.numeric(tmp3[1,2] - tmp3[2,2]) 
      res = rbind(res,c(j,FC,pvalue,i))
    }
  }
  
  res = as.data.frame(res)
  colnames(res) = c("Immune","FC","pvalue","subtype")
  res$FC = as.numeric(res$FC)
  res$pvalue = as.numeric(res$pvalue)
  res$Significance = ifelse(res$pvalue <= 0.05,"P.value <= 0.05","P.value > 0.05")
  
  p = ggplot(data = res)+
    ggtitle("Immune-related signatures (Protein)")+
    theme_bw()+
    geom_point(aes(x = Immune,y = subtype,size = -log10(pvalue),fill = FC,color = Significance),shape = 21,)+
    scale_color_manual(values = c("P.value <= 0.05" = "red","P.value > 0.05" = "grey"))+
    # scale_fill_viridis(limits = c(-max(res$FC),max(res$FC)),breaks = c(-max(res$FC),max(res$FC)),labels = c(paste(genes,"Low ↑↑"),paste(genes,"High ↑↑")))+
    scale_fill_gradient2(low = "blue3",mid = "white", high = "yellow3",limits = c(-max(abs(res$FC)),max(abs(res$FC))),breaks = c(-max(abs(res$FC)),max(abs(res$FC))),labels = c(paste(genes,"Low (Up)"),paste(genes,"High (Up)")))+
    theme(axis.text.x = element_text(angle = 30,size = 12,hjust = 1,face = "bold"),plot.title = element_text(size = 22,face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12,face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.title = element_blank(),
          axis.text.y = element_text(face = "bold",size = 15))
  
  return(p)
  
  
}

DEG_CPTAC_rna_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["rna"]]),patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    tmp_data = TCGA[[cancer_type]][["rna"]][,c(UP_name,DOWN_name)]

    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][["rna"]]),mut_patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    
    tmp_data = TCGA[[cancer_type]][["rna"]][,c(UP_name,DOWN_name)]
    
    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["rna"]]),wt_patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    tmp_data = TCGA[[cancer_type]][["rna"]][,c(UP_name,DOWN_name)]

    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }
}
DEG_CPTAC_rna_subtype = function(TCGA,min.pct,tmp_data){
  
  UP_name = tmp_data$UP
  DOWN_name = tmp_data$DOWN
  tmp_data = tmp_data$tmp_data
  
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),levels = c("High","Low"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("High","Low")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(High-Low,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
  
}

DEG_CPTAC_protein_subtype_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(colnames(TCGA[[cancer_type]][["protein"]]),patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id[tmp_mean_exp > cutoff]
    DOWN_name = pat_id[tmp_mean_exp <= cutoff]
    
    
    tmp_data = TCGA[[cancer_type]][["protein"]][,c(UP_name,DOWN_name)]
    
    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(colnames(TCGA[[cancer_type]][["protein"]]),mut_patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut])
    cutoff = median(tmp_mean_exp)
    UP_name = mut[tmp_mean_exp > cutoff]
    DOWN_name = mut[tmp_mean_exp <= cutoff]
    
    
    tmp_data = TCGA[[cancer_type]][["protein"]][,c(UP_name,DOWN_name)]

    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }else{
    
    wt = intersect(colnames(TCGA[[cancer_type]][["protein"]]),wt_patient_id)
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt])
    cutoff = median(tmp_mean_exp)
    UP_name = wt[tmp_mean_exp > cutoff]
    DOWN_name = wt[tmp_mean_exp <= cutoff]
    
    tmp_data = TCGA[[cancer_type]][["protein"]][,c(UP_name,DOWN_name)]

    return(list("tmp_data"=tmp_data,"UP"=UP_name,"DOWN"=DOWN_name))
    
  }
}
DEG_CPTAC_protein_subtype = function(TCGA,min.pct,tmp_data){
  
  UP_name = tmp_data$UP
  DOWN_name = tmp_data$DOWN
  tmp_data = tmp_data$tmp_data
  
  m = apply(tmp_data,1,function(x){sum(x>0)>ncol(tmp_data)*min.pct})
  tmp_data = tmp_data[m,]
  group = factor(c(rep("High",length(UP_name)),rep("Low",length(DOWN_name))),levels = c("High","Low"))
  design = model.matrix(~0+group)
  rownames(design) = colnames(tmp_data)
  colnames(design) = c("High","Low")
  
  fit = lmFit(tmp_data,design)
  contr = makeContrasts(High-Low,levels = design)
  diff = contrasts.fit(fit,contr)
  diff = eBayes(diff)
  tab = topTable(diff, sort.by = "P", n = Inf)
  return(tab)
 
}

CPTAC_survival_subtype_rna_fun_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,patient_id)
    pat_id2 = intersect(pat_id,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% pat_id2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id2])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id2[tmp_mean_exp > cutoff]
    DOWN_name = pat_id2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[tmp_data$PatientID %in% UP_name,"groups"] = "High"
    tmp_data[tmp_data$PatientID %in% DOWN_name,"groups"] = "low"
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,mut_patient_id)
    mut2 = intersect(mut,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% mut2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut2])
    cutoff = median(tmp_mean_exp)
    UP_name = mut2[tmp_mean_exp > cutoff]
    DOWN_name = mut2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[tmp_data$PatientID %in% UP_name,"groups"] = "High"
    tmp_data[tmp_data$PatientID %in% DOWN_name,"groups"] = "low"
    
  }else{
    
    wt = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,wt_patient_id)
    wt2 = intersect(wt,colnames(TCGA[[cancer_type]][["rna"]]))
    tmp_data = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% wt2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt2])
    cutoff = median(tmp_mean_exp)
    UP_name = wt2[tmp_mean_exp > cutoff]
    DOWN_name = wt2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[tmp_data$PatientID %in% UP_name,"groups"] = "High"
    tmp_data[tmp_data$PatientID %in% DOWN_name,"groups"] = "low"
  }
  
  tmp_data = tmp_data[,-1]
  return(tmp_data)
}
CPTAC_survival_subtype_rna_fun = function(TCGA,tmp_data){
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title="Genes expression(RNA)",
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  p = (p1$plot/p1$cumevents)+plot_layout(design = layout)
  return(p)
  
}

CPTAC_survival_subtype_protein_fun_tmpdata = function(TCGA,cancer_type,genes,MT_OR_WT,mut_patient_id,wt_patient_id){
  if(MT_OR_WT == "All"){
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,patient_id)
    pat_id2 = intersect(pat_id,colnames(TCGA[[cancer_type]][["protein"]]))
    tmp_data = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% pat_id2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id2])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id2[tmp_mean_exp > cutoff]
    DOWN_name = pat_id2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[tmp_data$PatientID %in% UP_name,"groups"] = "High"
    tmp_data[tmp_data$PatientID %in% DOWN_name,"groups"] = "low"
    
  }else if(MT_OR_WT == "Mutation"){
    
    mut = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,mut_patient_id)
    mut2 = intersect(mut,colnames(TCGA[[cancer_type]][["protein"]]))
    tmp_data = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% mut2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut2])
    cutoff = median(tmp_mean_exp)
    UP_name = mut2[tmp_mean_exp > cutoff]
    DOWN_name = mut2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[tmp_data$PatientID %in% UP_name,"groups"] = "High"
    tmp_data[tmp_data$PatientID %in% DOWN_name,"groups"] = "low"
    
  }else{
    
    wt = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,wt_patient_id)
    wt2 = intersect(wt,colnames(TCGA[[cancer_type]][["protein"]]))
    tmp_data = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% wt2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt2])
    cutoff = median(tmp_mean_exp)
    UP_name = wt2[tmp_mean_exp > cutoff]
    DOWN_name = wt2[tmp_mean_exp <= cutoff]
    
    tmp_data$groups = NA
    tmp_data[tmp_data$PatientID %in% UP_name,"groups"] = "High"
    tmp_data[tmp_data$PatientID %in% DOWN_name,"groups"] = "low"
  }
  
  tmp_data = tmp_data[,-1]
  return(tmp_data)
}
CPTAC_survival_subtype_protein_fun = function(TCGA,tmp_data){
  
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp_data[,"groups"],data = tmp_data))
  p1 = ggsurvplot(fit,data = tmp_data,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)",
                  legend.title="Genes expression(protein)",
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20)))
  
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  p = (p1$plot/p1$cumevents)+plot_layout(design = layout)
  return(p)
  
}

CPTAC_survival_subtype_fun_tmpdata_compare = function(TCGA,cancer_type,genes,mut_patient_id,wt_patient_id){
  
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,patient_id)
    pat_id2 = intersect(pat_id,colnames(TCGA[[cancer_type]][["rna"]]))
    All = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% pat_id2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,pat_id2])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id2[tmp_mean_exp > cutoff]
    DOWN_name = pat_id2[tmp_mean_exp <= cutoff]
    
    All$groups = NA
    All[All$PatientID %in% UP_name,"groups"] = "High"
    All[All$PatientID %in% DOWN_name,"groups"] = "low"
    All$omics = "RNA"
    All$module = "All"
    
    
    mut = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,mut_patient_id)
    mut2 = intersect(mut,colnames(TCGA[[cancer_type]][["rna"]]))
    Mutation = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% mut2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,mut2])
    cutoff = median(tmp_mean_exp)
    UP_name = mut2[tmp_mean_exp > cutoff]
    DOWN_name = mut2[tmp_mean_exp <= cutoff]
    
    Mutation$groups = NA
    Mutation[Mutation$PatientID %in% UP_name,"groups"] = "High"
    Mutation[Mutation$PatientID %in% DOWN_name,"groups"] = "low"
    Mutation$omics = "RNA"
    Mutation$module = paste(genes,"Mutation")

    
    wt = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,wt_patient_id)
    wt2 = intersect(wt,colnames(TCGA[[cancer_type]][["rna"]]))
    Wildtype = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% wt2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$rna[genes,wt2])
    cutoff = median(tmp_mean_exp)
    UP_name = wt2[tmp_mean_exp > cutoff]
    DOWN_name = wt2[tmp_mean_exp <= cutoff]
    
    Wildtype$groups = NA
    Wildtype[Wildtype$PatientID %in% UP_name,"groups"] = "High"
    Wildtype[Wildtype$PatientID %in% DOWN_name,"groups"] = "low"
    Wildtype$omics = "RNA"
    Wildtype$module = paste(genes,"Wildtype")
  
    RNA = rbind(All,Mutation,Wildtype)[,-1]
    
    patient_id = unique(TCGA[[cancer_type]][["maf"]]@data$TCGA)
    pat_id = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,patient_id)
    pat_id2 = intersect(pat_id,colnames(TCGA[[cancer_type]][["protein"]]))
    All = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% pat_id2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,pat_id2])
    cutoff = median(tmp_mean_exp)
    UP_name = pat_id2[tmp_mean_exp > cutoff]
    DOWN_name = pat_id2[tmp_mean_exp <= cutoff]
    
    All$groups = NA
    All[All$PatientID %in% UP_name,"groups"] = "High"
    All[All$PatientID %in% DOWN_name,"groups"] = "low"
    All$omics = "Protein"
    All$module = "All"
    
    mut = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,mut_patient_id)
    mut2 = intersect(mut,colnames(TCGA[[cancer_type]][["protein"]]))
    Mutation = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% mut2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,mut2])
    cutoff = median(tmp_mean_exp)
    UP_name = mut2[tmp_mean_exp > cutoff]
    DOWN_name = mut2[tmp_mean_exp <= cutoff]
    
    Mutation$groups = NA
    Mutation[Mutation$PatientID %in% UP_name,"groups"] = "High"
    Mutation[Mutation$PatientID %in% DOWN_name,"groups"] = "low"
    Mutation$omics = "Protein"
    Mutation$module = paste(genes,"Mutation")
    
    wt = intersect(TCGA[[cancer_type]][["Survival"]]$PatientID,wt_patient_id)
    wt2 = intersect(wt,colnames(TCGA[[cancer_type]][["protein"]]))
    Wildtype = TCGA[[cancer_type]][["Survival"]][TCGA[[cancer_type]][["Survival"]]$PatientID %in% wt2,]
    
    tmp_mean_exp = colMeans(TCGA[[cancer_type]]$protein[genes,wt2])
    cutoff = median(tmp_mean_exp)
    UP_name = wt2[tmp_mean_exp > cutoff]
    DOWN_name = wt2[tmp_mean_exp <= cutoff]
    
    Wildtype$groups = NA
    Wildtype[Wildtype$PatientID %in% UP_name,"groups"] = "High"
    Wildtype[Wildtype$PatientID %in% DOWN_name,"groups"] = "low"
    Wildtype$omics = "Protein"
    Wildtype$module = paste(genes,"Wildtype")
  
    Protein = rbind(All,Mutation,Wildtype)[,-1]
  
  return(rbind(RNA,Protein))
  
}
CPTAC_survival_subtype_fun_compare = function(tmp_data_compare,genes){
  
  layout = "
      AAAA
      AAAA
      AAAA
      AAAA
      BBBB
      "
  
  tmp2 = tmp_data_compare[ tmp_data_compare$omics == "RNA" & tmp_data_compare$module == "All",]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p1 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)\n(All)",
                  legend.title=paste(genes,"(RNA)"),
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.title = element_text(size = 15,face = "bold")))
  
  tmp2 = tmp_data_compare[ tmp_data_compare$omics == "RNA" & grepl(pattern = "Mutation",tmp_data_compare$module),]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p2 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title=paste("Overall Survival (OS)","\n(",tmp2$module[1],")",sep = ""),
                  legend.title=paste(genes,"(RNA)"),
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.title = element_text(size = 15,face = "bold")))
  
  tmp2 = tmp_data_compare[ tmp_data_compare$omics == "RNA" & grepl(pattern = "Wildtype",tmp_data_compare$module),]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p3 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title=paste("Overall Survival (OS)","\n(",tmp2$module[1],")",sep = ""),
                  legend.title=paste(genes,"(RNA)"),
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.title = element_text(size = 15,face = "bold")))
  

  RNA = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout) | (p3$plot/p3$cumevents)+plot_layout(design = layout)
  
  tmp2 = tmp_data_compare[ tmp_data_compare$omics == "Protein" & tmp_data_compare$module == "All",]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p1 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title="Overall Survival (OS)\n(All)",
                  legend.title=paste(genes,"(Protein)"),
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.title = element_text(size = 15,face = "bold")))
  
  tmp2 = tmp_data_compare[ tmp_data_compare$omics == "Protein" & grepl(pattern = "Mutation",tmp_data_compare$module),]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p2 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title=paste("Overall Survival (OS)","\n(",tmp2$module[1],")",sep = ""),
                  legend.title=paste(genes,"(Protein)"),
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.title = element_text(size = 15,face = "bold")))
  
  tmp2 = tmp_data_compare[ tmp_data_compare$omics == "Protein" & grepl(pattern = "Wildtype",tmp_data_compare$module),]
  fit <- do.call(survfit, list(Surv(OS.time,OS)~tmp2[,"groups"],data = tmp2))
  p3 = ggsurvplot(fit,data = tmp2,
                  xlab = "day",
                  pval = T,
                  risk.table = F,
                  cumevents = T,
                  palette = c("red2", "blue4"),
                  title=paste("Overall Survival (OS)","\n(",tmp2$module[1],")",sep = ""),
                  legend.title=paste(genes,"(Protein)"),
                  legend.labs = c("High","low"),
                  ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5,size = 20),legend.title = element_text(size = 15,face = "bold")))
  
  
  Protein = (p1$plot/p1$cumevents)+plot_layout(design = layout) | (p2$plot/p2$cumevents)+plot_layout(design = layout) | (p3$plot/p3$cumevents)+plot_layout(design = layout)
  
  p = RNA / Protein
  return(p)
}
################################# Search gene ###################################
OS_single_forest = function(gene){
  
  cohort = vector()
  wildtype = vector()
  mutation = vector()
  HR_mean = vector()
  HR_lower = vector()
  HR_upper = vector()
  HR_text = vector()
  logtest = vector()
  waldtest = vector()
  for(i in setdiff(names(datasets),c("dataset1","dataset3","dataset6"))){
    
    if(rowSums(is.na(ref_total_OS_single[[i]][gene,])) == 0){
      
      cohort = c(cohort,dataset_name2[i])
      wildtype = c(wildtype,ref_total_OS_single[[i]][gene,"Patient(Wildtype)"])
      mutation = c(mutation,ref_total_OS_single[[i]][gene,"Patient(Mutation)"])
      
      HR_mean = c(HR_mean,round(ref_total_OS_single[[i]][gene,"HR(OS)"],digits = 3))
      HR_lower = c(HR_lower,round(ref_total_OS_single[[i]][gene,"Lower(95%, OS)"],digits = 3))
      HR_upper = c(HR_upper,round(ref_total_OS_single[[i]][gene,"Upper(95%, OS)"],digits = 3))
      HR_text = c(HR_text,paste(round(ref_total_OS_single[[i]][gene,"HR(OS)"],digits = 3),"(",round(ref_total_OS_single[[i]][gene,"Lower(95%, OS)"],digits = 3),"-",round(ref_total_OS_single[[i]][gene,"Upper(95%, OS)"],digits = 3),")"))
      logtest = c(logtest,round(ref_total_OS_single[[i]][gene,"Log_rank_test(OS)"],digits = 3))
      waldtest = c(waldtest,round(ref_total_OS_single[[i]][gene,"Wald_test(OS)"],digits = 3))
      
      
      
      
    }
    
  }
  
  ord = order(HR_mean)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"Wild type",wildtype[ord],sum(wildtype),sum(wildtype)),
    c(NA,"Mutation",mutation[ord],sum(mutation),sum(mutation)),
    c(NA,"HR(95%CI)",HR_text[ord],NA,NA)
    #                   c(NA,"P value(Log rank)",logtest[ord],NA),
    #                   c(NA,"P value(wald text)",waldtest[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,HR_mean[ord],NA,NA),"lower" = c(NA,NA,HR_lower[ord],NA,NA),"upper" = c(NA,NA,HR_upper[ord],NA,NA))
  # cochrane_from_rmeta[ cochrane_from_rmeta$mean %in% c(0,Inf),] = NA
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "HR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>% 
    forestplot(labeltext = tabletext2, clip = c(minlow-minlow/10,maxup+maxup/10),
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("Overall Survival(OS)","\n",gene),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE, 
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}

PFS_single_forest = function(gene){
  
  cohort = vector()
  wildtype = vector()
  mutation = vector()
  HR_mean = vector()
  HR_lower = vector()
  HR_upper = vector()
  HR_text = vector()
  logtest = vector()
  waldtest = vector()
  for(i in setdiff(names(datasets),c("dataset1","dataset3","dataset6"))){
    if(rowSums(is.na(ref_total_PFS_single[[i]][gene,])) == 0){
      cohort = c(cohort,dataset_name2[i])
      wildtype = c(wildtype,ref_total_PFS_single[[i]][gene,"Patient(Wildtype)"])
      mutation = c(mutation,ref_total_PFS_single[[i]][gene,"Patient(Mutation)"])
      
      HR_mean = c(HR_mean,round(ref_total_PFS_single[[i]][gene,"HR(PFS)"],digits = 3))
      HR_lower = c(HR_lower,round(ref_total_PFS_single[[i]][gene,"Lower(95%, PFS)"],digits = 3))
      HR_upper = c(HR_upper,round(ref_total_PFS_single[[i]][gene,"Upper(95%, PFS)"],digits = 3))
      HR_text = c(HR_text,paste(round(ref_total_PFS_single[[i]][gene,"HR(PFS)"],digits = 3),"(",round(ref_total_PFS_single[[i]][gene,"Lower(95%, PFS)"],digits = 3),"-",round(ref_total_PFS_single[[i]][gene,"Upper(95%, PFS)"],digits = 3),")"))
      logtest = c(logtest,round(ref_total_PFS_single[[i]][gene,"Log_rank_test(PFS)"],digits = 3))
      waldtest = c(waldtest,round(ref_total_PFS_single[[i]][gene,"Wald_test(PFS)"],digits = 3))
    }
    
    
    
    
    
    
  }
  ord = order(HR_mean)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"Wild type",wildtype[ord],sum(wildtype),sum(wildtype)),
    c(NA,"Mutation",mutation[ord],sum(mutation),sum(mutation)),
    c(NA,"HR(95%CI)",HR_text[ord],NA,NA)
    #                   c(NA,"P value(Log rank)",logtest[ord],NA),
    #                   c(NA,"P value(wald text)",waldtest[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,HR_mean[ord],NA,NA),"lower" = c(NA,NA,HR_lower[ord],NA,NA),"upper" = c(NA,NA,HR_upper[ord],NA,NA))
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "HR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>% 
    forestplot(labeltext = tabletext2, clip = c(minlow-minlow/10,maxup+maxup/10),
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("Progress Free Survival(PFS)","\n",gene),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE, 
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  return(p)
}

RECIST_single_forest = function(gene){
  
  cohort = vector()
  SDPD = vector()
  CRPR = vector()
  OR_mean = vector()
  OR_lower = vector()
  OR_upper = vector()
  OR_text = vector()
  midp.exact = vector()
  chi.square = vector()
  fisher.exact = vector()
  for(i in setdiff(names(datasets),c("dataset1","dataset3","dataset6"))){
    
    if(rowSums(is.na(ref_total_RECIST_single[[i]][gene,])) == 0){
      
      cohort = c(cohort,dataset_name2[i])
      SDPD = c(SDPD,ref_total_RECIST_single[[i]][gene,"SD & PD(Mutation/total)"])
      CRPR = c(CRPR,ref_total_RECIST_single[[i]][gene,"CR & PR(Mutation/total)"])
      
      OR_mean = c(OR_mean,ref_total_RECIST_single[[i]][gene,"OR"])
      OR_lower = c(OR_lower,ref_total_RECIST_single[[i]][gene,"Lower(95%)"])
      OR_upper = c(OR_upper,ref_total_RECIST_single[[i]][gene,"Upper(95%)"])
      OR_text = c(OR_text,paste(round(ref_total_RECIST_single[[i]][gene,"OR"],digits = 3),"(",round(ref_total_RECIST_single[[i]][gene,"Lower(95%)"],digits = 3),"-",round(ref_total_RECIST_single[[i]][gene,"Upper(95%)"],digits = 3),")"))
      midp.exact = c(midp.exact,round(ref_total_RECIST_single[[i]][gene,"midp.exact"],digits = 3))
      chi.square = c(chi.square,round(ref_total_RECIST_single[[i]][gene,"chi.square"],digits = 3))
      fisher.exact = c(fisher.exact,round(ref_total_RECIST_single[[i]][gene,"fisher.exact"],digits = 3))
      
      
      
      
    }
    
  }
  
  
  ord = order(OR_mean,decreasing = T)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"SD & PD\n(Mutation/total)",SDPD[ord],NA,NA),
    c(NA,"CR & PR\n(Mutation/total)",CRPR[ord],NA,NA),
    c(NA,"OR(95%CI)",OR_text[ord],NA,NA)
    #     c(NA,"P value\n(midp.exact)",midp.exact[ord],NA,NA),
    #     c(NA,"P value\n(chi.square)",chi.square[ord],NA,NA),
    #     c(NA,"P value\n(fisher.exact)",fisher.exact[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,OR_mean[ord],NA,NA),"lower" = c(NA,NA,OR_lower[ord],NA,NA),"upper" = c(NA,NA,OR_upper[ord],NA,NA))
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "OR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>%
    forestplot(labeltext = tabletext2,clip = c(minlow-minlow/10,maxup+maxup/10),
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("RECIST","\n",gene),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE,
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}

RESPONSE_single_forest = function(gene){
  
  cohort = vector()
  NCB = vector()
  DCB = vector()
  OR_mean = vector()
  OR_lower = vector()
  OR_upper = vector()
  OR_text = vector()
  midp.exact = vector()
  chi.square = vector()
  fisher.exact = vector()
  for(i in setdiff(names(datasets),c("dataset1","dataset3","dataset6"))){
    
    if(rowSums(is.na(ref_total_RESPONSE_single[[i]][gene,])) == 0){
      
      cohort = c(cohort,dataset_name2[i])
      NCB = c(NCB,ref_total_RESPONSE_single[[i]][gene,"NCB(Mutation/total)"])
      DCB = c(DCB,ref_total_RESPONSE_single[[i]][gene,"DCB(Mutation/total)"])
      
      OR_mean = c(OR_mean,ref_total_RESPONSE_single[[i]][gene,"OR"])
      OR_lower = c(OR_lower,ref_total_RESPONSE_single[[i]][gene,"Lower(95%)"])
      OR_upper = c(OR_upper,ref_total_RESPONSE_single[[i]][gene,"Upper(95%)"])
      OR_text = c(OR_text,paste(round(ref_total_RESPONSE_single[[i]][gene,"OR"],digits = 3),"(",round(ref_total_RESPONSE_single[[i]][gene,"Lower(95%)"],digits = 3),"-",round(ref_total_RESPONSE_single[[i]][gene,"Upper(95%)"],digits = 3),")"))
      midp.exact = c(midp.exact,round(ref_total_RESPONSE_single[[i]][gene,"midp.exact"],digits = 3))
      chi.square = c(chi.square,round(ref_total_RESPONSE_single[[i]][gene,"chi.square"],digits = 3))
      fisher.exact = c(fisher.exact,round(ref_total_RESPONSE_single[[i]][gene,"fisher.exact"],digits = 3))
      
      
      
      
    }
    
  }
  ord = order(OR_mean,decreasing = T)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"NCB\n(Mutation/total)",NCB[ord],NA,NA),
    c(NA,"DCB\n(Mutation/total)",DCB[ord],NA,NA),
    c(NA,"OR(95%CI)",OR_text[ord],NA,NA)
    #     c(NA,"P value\n(midp.exact)",midp.exact[ord],NA,NA),
    #     c(NA,"P value\n(chi.square)",chi.square[ord],NA,NA),
    #     c(NA,"P value\n(fisher.exact)",fisher.exact[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,OR_mean[ord],NA,NA),"lower" = c(NA,NA,OR_lower[ord],NA,NA),"upper" = c(NA,NA,OR_upper[ord],NA,NA))
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "OR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>%
    forestplot(labeltext = tabletext2,clip = c(minlow-minlow/10,maxup+maxup/10),
               
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("Clinical Benefit","\n",gene),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE,
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}

ref_single_Immune_Infiltration_heatmap = function(gene){

  All_diff = vector()
  All_pvalue = vector()
  for(i in setdiff(names(ref_total_immune_infiltrating_single),c("dataset1","dataset3","dataset6"))){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in inf_names){

      diff = ref_total_immune_infiltrating_single[[i]][[j]][gene,"Score_median(Mutation)"] - ref_total_immune_infiltrating_single[[i]][[j]][gene,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,ref_total_immune_infiltrating_single[[i]][[j]][gene,"P_value"])

    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = dataset_name2[setdiff(names(ref_total_immune_infiltrating_single),c("dataset1","dataset3","dataset6"))]
  colnames(All_pvalue) = dataset_name2[setdiff(names(ref_total_immune_infiltrating_single),c("dataset1","dataset3","dataset6"))]

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 30,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

ref_single_Immune_pathway_heatmap = function(gene){

  All_diff = vector()
  All_pvalue = vector()
  for(i in setdiff(names(ref_total_immune_pathway_single),c("dataset1","dataset3","dataset6"))){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in sig_names){

      diff = ref_total_immune_pathway_single[[i]][[j]][gene,"Score_median(Mutation)"] - ref_total_immune_pathway_single[[i]][[j]][gene,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,ref_total_immune_pathway_single[[i]][[j]][gene,"P_value"])

    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = dataset_name2[setdiff(names(ref_total_immune_pathway_single),c("dataset1","dataset3","dataset6"))]
  colnames(All_pvalue) = dataset_name2[setdiff(names(ref_total_immune_pathway_single),c("dataset1","dataset3","dataset6"))]

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 30,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

TCGA_single_Immune_Infiltration_heatmap = function(gene){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(TCGA)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in inf_names){

      diff = TCGA_total_immune_infiltration_single[[i]][[j]][gene,"Score_median(Mutation)"] - TCGA_total_immune_infiltration_single[[i]][[j]][gene,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,TCGA_total_immune_infiltration_single[[i]][[j]][gene,"P_value"])

    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(TCGA)
  colnames(All_pvalue) = names(TCGA)

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

TCGA_single_Immune_pathway_heatmap = function(gene){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(TCGA)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in sig_names){

      diff = TCGA_total_immune_pathway_single[[i]][[j]][gene,"Score_median(Mutation)"] - TCGA_total_immune_pathway_single[[i]][[j]][gene,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,TCGA_total_immune_pathway_single[[i]][[j]][gene,"P_value"])

    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(TCGA)
  colnames(All_pvalue) = names(TCGA)

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

CPTAC_rna_single_Immune_Infiltration_heatmap = function(gene){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_inf_names = names(CPTAC_total_immune_infiltration_rna_single[[i]])
    for(j in inf_names){
      if(j %in% tmp_inf_names){
        diff = CPTAC_total_immune_infiltration_rna_single[[i]][[j]][gene,"Score_median(Mutation)"] - CPTAC_total_immune_infiltration_rna_single[[i]][[j]][gene,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_infiltration_rna_single[[i]][[j]][gene,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)

      }


    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(inf_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(inf_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}
CPTAC_protein_single_Immune_Infiltration_heatmap = function(gene){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_inf_names = names(CPTAC_total_immune_infiltration_protein_single[[i]])
    for(j in inf_names){
      if(j %in% tmp_inf_names){
        diff = CPTAC_total_immune_infiltration_protein_single[[i]][[j]][gene,"Score_median(Mutation)"] - CPTAC_total_immune_infiltration_protein_single[[i]][[j]][gene,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_infiltration_protein_single[[i]][[j]][gene,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)
      }


    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(inf_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(inf_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

CPTAC_rna_single_Immune_pathway_heatmap = function(gene){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_sig_names = names(CPTAC_total_immune_pathway_rna_single[[i]])
    for(j in sig_names){
      if(j %in% tmp_sig_names){
        diff = CPTAC_total_immune_pathway_rna_single[[i]][[j]][gene,"Score_median(Mutation)"] - CPTAC_total_immune_pathway_rna_single[[i]][[j]][gene,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_pathway_rna_single[[i]][[j]][gene,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)
      }


    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(sig_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(sig_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l,na.rm = T)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l,na.rm = T)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l,na.rm = T)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}
CPTAC_protein_single_Immune_pathway_heatmap = function(gene){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_sig_names = names(CPTAC_total_immune_pathway_protein_single[[i]])
    for(j in sig_names){
      if(j %in% tmp_sig_names){
        diff = CPTAC_total_immune_pathway_protein_single[[i]][[j]][gene,"Score_median(Mutation)"] - CPTAC_total_immune_pathway_protein_single[[i]][[j]][gene,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_pathway_protein_single[[i]][[j]][gene,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)
      }


    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(sig_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(sig_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l,na.rm = T)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l,na.rm = T)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l,na.rm = T)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)}

Gene_detial_information = function(gene){
  myheaders <- c("user-agent"="Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Mobile Safari/537.36 Edg/109.0.1518.78")
  response = GET(paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene,sep = ""),add_headers(.headers = myheaders))

  gene_title = read_html(response)%>% html_nodes(xpath="//section[@id='summaries']//h3")%>% html_text()
  gene_text = read_html(response)%>% html_nodes(xpath="//section[@id='summaries']//p")%>% html_text()

  len = length(gene_text)
  gene_title = gene_title[1:len]

  res = character()
  for( i in 1:len){
    gene_text[i] = gsub("\r|\n|  ","",gene_text[i])
    gene_title[i] = gsub("\r|\n|  ","",gene_title[i])
    tmp = paste("<h3>",gene_title[i],"</h3>","<p>",gene_text[i],"</p>",sep = "")
    res = paste(res,tmp,sep = "")
  }
  res
}


################################# Search pathway ################################
OS_pm_forest = function(pathway){
  
  cohort = vector()
  wildtype = vector()
  mutation = vector()
  HR_mean = vector()
  HR_lower = vector()
  HR_upper = vector()
  HR_text = vector()
  logtest = vector()
  waldtest = vector()
  for(i in setdiff(rownames(datasets_overview)[datasets_overview$Sequencing == "WES"],c("dataset1","dataset3","dataset6"))){
    
    if(rowSums(is.na(ref_total_OS_pm[[i]][pathway,])) == 0){
      
      cohort = c(cohort,dataset_name2[i])
      wildtype = c(wildtype,ref_total_OS_pm[[i]][pathway,"Patient(Wildtype)"])
      mutation = c(mutation,ref_total_OS_pm[[i]][pathway,"Patient(Mutation)"])
      
      HR_mean = c(HR_mean,round(ref_total_OS_pm[[i]][pathway,"HR(OS)"],digits = 3))
      HR_lower = c(HR_lower,round(ref_total_OS_pm[[i]][pathway,"Lower(95%, OS)"],digits = 3))
      HR_upper = c(HR_upper,round(ref_total_OS_pm[[i]][pathway,"Upper(95%, OS)"],digits = 3))
      HR_text = c(HR_text,paste(round(ref_total_OS_pm[[i]][pathway,"HR(OS)"],digits = 3),"(",round(ref_total_OS_pm[[i]][pathway,"Lower(95%, OS)"],digits = 3),"-",round(ref_total_OS_pm[[i]][pathway,"Upper(95%, OS)"],digits = 3),")"))
      logtest = c(logtest,round(ref_total_OS_pm[[i]][pathway,"Log_rank_test(OS)"],digits = 3))
      waldtest = c(waldtest,round(ref_total_OS_pm[[i]][pathway,"Wald_test(OS)"],digits = 3))
      
      
      
      
    }
    
  }
  ord = order(HR_mean)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"Wild type",wildtype[ord],sum(wildtype),sum(wildtype)),
    c(NA,"Mutation",mutation[ord],sum(mutation),sum(mutation)),
    c(NA,"HR(95%CI)",HR_text[ord],NA,NA)
    #                   c(NA,"P value(Log rank)",logtest[ord],NA),
    #                   c(NA,"P value(wald text)",waldtest[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,HR_mean[ord],NA,NA),"lower" = c(NA,NA,HR_lower[ord],NA,NA),"upper" = c(NA,NA,HR_upper[ord],NA,NA))
  # cochrane_from_rmeta[ cochrane_from_rmeta$mean %in% c(0,Inf),] = NA
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "HR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>% 
    forestplot(labeltext = tabletext2, clip = c(minlow-minlow/10,maxup+maxup/10),
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("Overall Survival(OS)","\n",pathway),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE, 
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}
PFS_pm_forest = function(pathway){
  
  cohort = vector()
  wildtype = vector()
  mutation = vector()
  HR_mean = vector()
  HR_lower = vector()
  HR_upper = vector()
  HR_text = vector()
  logtest = vector()
  waldtest = vector()
  for(i in setdiff(rownames(datasets_overview)[datasets_overview$Sequencing == "WES"],c("dataset1","dataset3","dataset6"))){
    if(rowSums(is.na(ref_total_PFS_pm[[i]][pathway,])) == 0){
      cohort = c(cohort,dataset_name2[i])
      wildtype = c(wildtype,ref_total_PFS_pm[[i]][pathway,"Patient(Wildtype)"])
      mutation = c(mutation,ref_total_PFS_pm[[i]][pathway,"Patient(Mutation)"])
      
      HR_mean = c(HR_mean,round(ref_total_PFS_pm[[i]][pathway,"HR(PFS)"],digits = 3))
      HR_lower = c(HR_lower,round(ref_total_PFS_pm[[i]][pathway,"Lower(95%, PFS)"],digits = 3))
      HR_upper = c(HR_upper,round(ref_total_PFS_pm[[i]][pathway,"Upper(95%, PFS)"],digits = 3))
      HR_text = c(HR_text,paste(round(ref_total_PFS_pm[[i]][pathway,"HR(PFS)"],digits = 3),"(",round(ref_total_PFS_pm[[i]][pathway,"Lower(95%, PFS)"],digits = 3),"-",round(ref_total_PFS_pm[[i]][pathway,"Upper(95%, PFS)"],digits = 3),")"))
      logtest = c(logtest,round(ref_total_PFS_pm[[i]][pathway,"Log_rank_test(PFS)"],digits = 3))
      waldtest = c(waldtest,round(ref_total_PFS_pm[[i]][pathway,"Wald_test(PFS)"],digits = 3))
    }
    
    
    
    
    
    
  }
  ord = order(HR_mean)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"Wild type",wildtype[ord],sum(wildtype),sum(wildtype)),
    c(NA,"Mutation",mutation[ord],sum(mutation),sum(mutation)),
    c(NA,"HR(95%CI)",HR_text[ord],NA,NA)
    #                   c(NA,"P value(Log rank)",logtest[ord],NA),
    #                   c(NA,"P value(wald text)",waldtest[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,HR_mean[ord],NA,NA),"lower" = c(NA,NA,HR_lower[ord],NA,NA),"upper" = c(NA,NA,HR_upper[ord],NA,NA))
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "HR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>% 
    forestplot(labeltext = tabletext2, clip = c(minlow-minlow/10,maxup+maxup/10),
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("Progress Free Survival(PFS)","\n",pathway),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE, 
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}
RECIST_pm_forest = function(pathway){
  
  cohort = vector()
  SDPD = vector()
  CRPR = vector()
  OR_mean = vector()
  OR_lower = vector()
  OR_upper = vector()
  OR_text = vector()
  midp.exact = vector()
  chi.square = vector()
  fisher.exact = vector()
  for(i in setdiff(rownames(datasets_overview)[datasets_overview$Sequencing == "WES"],c("dataset1","dataset3","dataset6"))){
    
    if(rowSums(is.na(ref_total_RECIST_pm[[i]][pathway,])) == 0){
      
      cohort = c(cohort,dataset_name2[i])
      SDPD = c(SDPD,ref_total_RECIST_pm[[i]][pathway,"SD & PD(Mutation/total)"])
      CRPR = c(CRPR,ref_total_RECIST_pm[[i]][pathway,"CR & PR(Mutation/total)"])
      
      OR_mean = c(OR_mean,ref_total_RECIST_pm[[i]][pathway,"OR"])
      OR_lower = c(OR_lower,ref_total_RECIST_pm[[i]][pathway,"Lower(95%)"])
      OR_upper = c(OR_upper,ref_total_RECIST_pm[[i]][pathway,"Upper(95%)"])
      OR_text = c(OR_text,paste(round(ref_total_RECIST_pm[[i]][pathway,"OR"],digits = 3),"(",round(ref_total_RECIST_pm[[i]][pathway,"Lower(95%)"],digits = 3),"-",round(ref_total_RECIST_pm[[i]][pathway,"Upper(95%)"],digits = 3),")"))
      midp.exact = c(midp.exact,round(ref_total_RECIST_pm[[i]][pathway,"midp.exact"],digits = 3))
      chi.square = c(chi.square,round(ref_total_RECIST_pm[[i]][pathway,"chi.square"],digits = 3))
      fisher.exact = c(fisher.exact,round(ref_total_RECIST_pm[[i]][pathway,"fisher.exact"],digits = 3))
      
      
      
      
    }
    
  }
  ord = order(OR_mean,decreasing = T)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"SD & PD\n(Mutation/total)",SDPD[ord],NA,NA),
    c(NA,"CR & PR\n(Mutation/total)",CRPR[ord],NA,NA),
    c(NA,"OR(95%CI)",OR_text[ord],NA,NA)
    #     c(NA,"P value\n(midp.exact)",midp.exact[ord],NA,NA),
    #     c(NA,"P value\n(chi.square)",chi.square[ord],NA,NA),
    #     c(NA,"P value\n(fisher.exact)",fisher.exact[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,OR_mean[ord],NA,NA),"lower" = c(NA,NA,OR_lower[ord],NA,NA),"upper" = c(NA,NA,OR_upper[ord],NA,NA))
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "OR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>%
    forestplot(labeltext = tabletext2,clip = c(minlow-minlow/10,maxup+maxup/10),
               
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("RECIST","\n",pathway),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE,
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}
RESPONSE_pm_forest = function(pathway){
  
  cohort = vector()
  NCB = vector()
  DCB = vector()
  OR_mean = vector()
  OR_lower = vector()
  OR_upper = vector()
  OR_text = vector()
  midp.exact = vector()
  chi.square = vector()
  fisher.exact = vector()
  for(i in setdiff(rownames(datasets_overview)[datasets_overview$Sequencing == "WES"],c("dataset1","dataset3","dataset6"))){
    
    if(rowSums(is.na(ref_total_RESPONSE_pm[[i]][pathway,])) == 0){
      
      cohort = c(cohort,dataset_name2[i])
      NCB = c(NCB,ref_total_RESPONSE_pm[[i]][pathway,"NCB(Mutation/total)"])
      DCB = c(DCB,ref_total_RESPONSE_pm[[i]][pathway,"DCB(Mutation/total)"])
      
      OR_mean = c(OR_mean,ref_total_RESPONSE_pm[[i]][pathway,"OR"])
      OR_lower = c(OR_lower,ref_total_RESPONSE_pm[[i]][pathway,"Lower(95%)"])
      OR_upper = c(OR_upper,ref_total_RESPONSE_pm[[i]][pathway,"Upper(95%)"])
      OR_text = c(OR_text,paste(round(ref_total_RESPONSE_pm[[i]][pathway,"OR"],digits = 3),"(",round(ref_total_RESPONSE_pm[[i]][pathway,"Lower(95%)"],digits = 3),"-",round(ref_total_RESPONSE_pm[[i]][pathway,"Upper(95%)"],digits = 3),")"))
      midp.exact = c(midp.exact,round(ref_total_RESPONSE_pm[[i]][pathway,"midp.exact"],digits = 3))
      chi.square = c(chi.square,round(ref_total_RESPONSE_pm[[i]][pathway,"chi.square"],digits = 3))
      fisher.exact = c(fisher.exact,round(ref_total_RESPONSE_pm[[i]][pathway,"fisher.exact"],digits = 3))
      
      
      
      
    }
    
  }
  ord = order(OR_mean,decreasing = T)
  tabletext = cbind(
    c(NA,"Cohort",cohort[ord],"Common effect model","Random effect model"),
    c(NA,"NCB\n(Mutation/total)",NCB[ord],NA,NA),
    c(NA,"DCB\n(Mutation/total)",DCB[ord],NA,NA),
    c(NA,"OR(95%CI)",OR_text[ord],NA,NA)
    #     c(NA,"P value\n(midp.exact)",midp.exact[ord],NA,NA),
    #     c(NA,"P value\n(chi.square)",chi.square[ord],NA,NA),
    #     c(NA,"P value\n(fisher.exact)",fisher.exact[ord],NA,NA)
  )
  cochrane_from_rmeta = data.frame("mean" = c(NA,NA,OR_mean[ord],NA,NA),"lower" = c(NA,NA,OR_lower[ord],NA,NA),"upper" = c(NA,NA,OR_upper[ord],NA,NA))
  
  tmp = cochrane_from_rmeta[3:(nrow(cochrane_from_rmeta)-2),]
  rownames(tmp) = tabletext[3:(nrow(cochrane_from_rmeta)-2),1]
  tmp$se = (log2(tmp$upper) - log2(tmp$lower))/(2*1.96)
  m = metagen(TE = log2(tmp$mean),seTE = tmp$se,studlab = rownames(tmp),sm = "OR", method.tau = "DL")
  k = summary(m)
  
  chtext = paste(round(exp(x = 1)^(k$common$TE),3),"[",round(exp(x = 1)^(k$common$lower),3),",",round(exp(x = 1)^(k$common$upper),3),"]")
  tabletext[nrow(tabletext)-1,4] = c(chtext)
  cochrane_from_rmeta[nrow(tabletext)-1,1:3] = c(round(exp(x = 1)^(k$common$TE),3),round(exp(x = 1)^(k$common$lower),3),round(exp(x = 1)^(k$common$upper),3))
  rhtext = paste(round(exp(x = 1)^(k$random$TE),3),"[",round(exp(x = 1)^(k$random$lower),3),",",round(exp(x = 1)^(k$random$upper),3),"]")
  tabletext[nrow(tabletext),4] = c(rhtext)
  cochrane_from_rmeta[nrow(tabletext),1:3] = c(round(exp(x = 1)^(k$random$TE),3),round(exp(x = 1)^(k$random$lower),3),round(exp(x = 1)^(k$random$upper),3))
  
  if(k$pval.Q < 0.1){
    tabletext2 = tabletext[-(nrow(tabletext)-1),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(random)",paste(round(k$w.random/sum(k$w.random),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.random,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.random,digits = 4),")")
    }
    
    
    
  }else{
    tabletext2 = tabletext[-(nrow(tabletext)),]
    cochrane_from_rmeta2 = cochrane_from_rmeta[-(nrow(tabletext)-1),]
    tabletext2 = cbind(tabletext2,c(NA,"Weight(common)",paste(round(k$w.common/sum(k$w.common),3)*100,"%",sep = ""),"100%"))
    
    if(round(k$pval.common,digits = 4)<0.0001){
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta < 0.0001)")
    }else{
      tabletext2[nrow(tabletext2),1] = paste(tabletext2[nrow(tabletext2),1],"(P-meta = ",round(k$pval.common,digits = 4),")")
    }
    
  }
  
  maxup = max(cochrane_from_rmeta2$upper[!is.infinite(cochrane_from_rmeta2$upper)],na.rm = T)
  minlow = min(cochrane_from_rmeta2$lower[cochrane_from_rmeta2$lower != 0],na.rm = T)
  
  tmp_list = list()
  tmp_list[["3"]] = gpar(col = "#000044")
  tmp_list[[as.character(nrow(cochrane_from_rmeta2) + 1)]] = gpar(col = "#000044")
  
  options(repr.plot.height = 10, repr.plot.width = 20)
  p = cochrane_from_rmeta2 %>%
    forestplot(labeltext = tabletext2,clip = c(minlow-minlow/10,maxup+maxup/10),
               
               lty.ci = 1,
               ci.vertices.height = 0.3,
               align="c",
               graphwidth = unit(8,"cm"),
               lwd.zero = gpar(lwd=2,lty=2),
               txt_gp = fpTxtGp(title = gpar(cex = 2),ticks=gpar(cex=1.1),summary=gpar(cex = 1.8),label = gpar(cex = 1.5)),
               boxsize = 0.5,
               graph.pos = 4,
               title = paste("Clinical Benefit","\n",pathway),
               is.summary = c(rep(TRUE, 2), rep(FALSE, nrow(cochrane_from_rmeta2)-3)),
               xlog = TRUE,
               hrzl_lines = tmp_list,
               vertices = TRUE,
               col = fpColors(box = "royalblue",
                              line = "darkblue",
                              summary = "royalblue")
    )
  print(p)
}
ref_pm_Immune_Infiltration_heatmap = function(pathway){

  All_diff = vector()
  All_pvalue = vector()
  tmp_dataset_names = setdiff(intersect(names(ref_total_immune_infiltrating_pm),
                                rownames(datasets_overview)[datasets_overview$Sequencing == "WES"]),
                              c("dataset1","dataset3","dataset6"))

  for(i in tmp_dataset_names){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in inf_names){

      diff = ref_total_immune_infiltrating_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - ref_total_immune_infiltrating_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,ref_total_immune_infiltrating_pm[[i]][[j]][pathway,"P_value"])

    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = dataset_name2[tmp_dataset_names]
  colnames(All_pvalue) = dataset_name2[tmp_dataset_names]

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 30,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}
ref_pm_Immune_pathway_heatmap = function(pathway){

  All_diff = vector()
  All_pvalue = vector()
  tmp_dataset_names = setdiff(intersect(names(ref_total_immune_pathway_pm),
                                rownames(datasets_overview)[datasets_overview$Sequencing == "WES"]),
                       c("dataset1","dataset3","dataset6"))
  for(i in tmp_dataset_names){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in sig_names){

      diff = ref_total_immune_pathway_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - ref_total_immune_pathway_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,ref_total_immune_pathway_pm[[i]][[j]][pathway,"P_value"])

    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = dataset_name2[tmp_dataset_names]
  colnames(All_pvalue) = dataset_name2[tmp_dataset_names]

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 30,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

TCGA_pm_Immune_Infiltration_heatmap = function(pathway){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(TCGA)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in inf_names){

      diff = TCGA_total_immune_infiltration_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - TCGA_total_immune_infiltration_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,TCGA_total_immune_infiltration_pm[[i]][[j]][pathway,"P_value"])

    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(TCGA)
  colnames(All_pvalue) = names(TCGA)

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}
TCGA_pm_Immune_pathway_heatmap = function(pathway){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(TCGA)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    for(j in sig_names){

      diff = TCGA_total_immune_pathway_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - TCGA_total_immune_pathway_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
      one_cancer_diff = c(one_cancer_diff,diff)
      one_cancer_pvalue = c(one_cancer_pvalue,TCGA_total_immune_pathway_pm[[i]][[j]][pathway,"P_value"])

    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(TCGA)
  colnames(All_pvalue) = names(TCGA)

  All_diff = All_diff[,!is.na(colSums(All_diff))]
  All_pvalue = All_pvalue[,!is.na(colSums(All_pvalue))]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

CPTAC_rna_pm_Immune_Infiltration_heatmap = function(pathway){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_inf_names = names(CPTAC_total_immune_infiltration_rna_pm[[i]])
    for(j in inf_names){
      if(j %in% tmp_inf_names){
        diff = CPTAC_total_immune_infiltration_rna_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - CPTAC_total_immune_infiltration_rna_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_infiltration_rna_pm[[i]][[j]][pathway,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)

      }


    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(inf_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(inf_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}
CPTAC_protein_pm_Immune_Infiltration_heatmap = function(pathway){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_inf_names = names(CPTAC_total_immune_infiltration_protein_pm[[i]])
    for(j in inf_names){
      if(j %in% tmp_inf_names){
        diff = CPTAC_total_immune_infiltration_protein_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - CPTAC_total_immune_infiltration_protein_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_infiltration_protein_pm[[i]][[j]][pathway,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)
      }


    }
    names(one_cancer_diff) = inf_names
    names(one_cancer_pvalue) = inf_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(inf_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(inf_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}

CPTAC_rna_pm_Immune_pathway_heatmap = function(pathway){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_sig_names = names(CPTAC_total_immune_pathway_rna_pm[[i]])
    for(j in sig_names){
      if(j %in% tmp_sig_names){
        diff = CPTAC_total_immune_pathway_rna_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - CPTAC_total_immune_pathway_rna_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_pathway_rna_pm[[i]][[j]][pathway,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)
      }


    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(sig_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(sig_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l,na.rm = T)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l,na.rm = T)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l,na.rm = T)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)
}
CPTAC_protein_pm_Immune_pathway_heatmap = function(pathway){
  All_diff = vector()
  All_pvalue = vector()
  for(i in names(CPTAC)){
    one_cancer_diff = vector()
    one_cancer_pvalue = vector()
    tmp_sig_names = names(CPTAC_total_immune_pathway_protein_pm[[i]])
    for(j in sig_names){
      if(j %in% tmp_sig_names){
        diff = CPTAC_total_immune_pathway_protein_pm[[i]][[j]][pathway,"Score_median(Mutation)"] - CPTAC_total_immune_pathway_protein_pm[[i]][[j]][pathway,"Score_median(Wildtype)"]
        one_cancer_diff = c(one_cancer_diff,diff)
        one_cancer_pvalue = c(one_cancer_pvalue,CPTAC_total_immune_pathway_protein_pm[[i]][[j]][pathway,"P_value"])
      }else{
        one_cancer_diff = c(one_cancer_diff,NA)
        one_cancer_pvalue = c(one_cancer_pvalue,NA)
      }


    }
    names(one_cancer_diff) = sig_names
    names(one_cancer_pvalue) = sig_names

    All_diff = cbind(All_diff,one_cancer_diff)
    All_pvalue = cbind(All_pvalue,one_cancer_pvalue)


  }

  colnames(All_diff) = names(CPTAC)
  colnames(All_pvalue) = names(CPTAC)

  All_diff = All_diff[,colSums(is.na(All_diff)) != length(sig_names)]
  All_pvalue = All_pvalue[,colSums(is.na(All_pvalue)) != length(sig_names)]


  p = Heatmap(All_diff,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",
              #             col = col_fun,
              cluster_rows = T,name = "DIFF",
              row_names_gp = gpar(fontface = "bold",cex = 1),
              column_names_gp = gpar(fontface = "bold",cex = 1),
              row_names_side = "left",column_names_rot = 45,column_names_centered = FALSE,
              layer_fun = function(j, i, x, y, w, h, fill) {
                v = pindex(All_pvalue, i, j)
                l = v <=0.05 & v > 0.01
                if(any(l,na.rm = T)){
                  grid.text("*",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v <= 0.01 & v > 0.001
                if(any(l,na.rm = T)){
                  grid.text("**",x[l], y[l],gp = gpar(fontsize = 20))
                }

                l = v < 0.001
                if(any(l,na.rm = T)){
                  grid.text("***",x[l], y[l],gp = gpar(fontsize = 20))
                }

              }
  )

  print(p)}
#######################################################################################################
volcano_plot = function(tab,FC,pvalue){
  tab$gene = rownames(tab)
  tab$color = NA
  tab$color[tab$logFC > FC & tab$P.Value < pvalue] = "UP"
  tab$color[tab$logFC < -FC & tab$P.Value < pvalue] = "DOWN"
  tab$color[is.na(tab$color)] = "NOT"
  tab$color = factor(tab$color,levels = c("DOWN","NOT","UP"))
  # p = ggplot(data = tab,aes(x = logFC,y = -log10(P.Value),color = color))+
  #     geom_point()+
  #     geom_vline(xintercept = c(-log2(FC),log2(FC)),linetype = "dashed",color = "gray")+
  #     geom_hline(yintercept = -log10(pvalue),linetype = "dashed",color = "gray")+
  #     scale_color_manual(values = c("blue","gray","red"))+
  #     theme_classic()
  # return(p)
  
  
  plotly::config(plot_ly(tab,x = ~logFC, y = ~-log10(P.Value),text = ~gene, 
                         type = 'scatter',  mode = 'markers',color = ~color,
                         colors = c("blue","gray","red")
                         ),
         displayModeBar = FALSE
         )
}

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
gseaScores <- getFromNamespace("gseaScores", "DOSE")

tableGrob2 <- function(d, p = NULL) {
  # has_package("gridExtra")
  d <- d[order(rownames(d)),]
  tp <- gridExtra::tableGrob(d,theme = gridExtra::ttheme_default(base_size = 20))
  if (is.null(p)) {
    return(tp)
  }
  
  # Fix bug: The 'group' order of lines and dots/path is different
  p_data <- ggplot_build(p)$data[[1]]
  # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  p_data <- p_data[order(p_data[["group"]]), ]
  pcol <- unique(p_data[["colour"]])
  ## This is fine too
  ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
  }
  return(tp)
}

my_gseaplot2 = function (x, geneSetID, title = "",self.Description = geneSetID, color = "green", base_size = 11, 
                         rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, title.size=22,
                         ES_geom = "line") 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), axis.title.y = element_text(size = 15,face = "bold"),
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
                                     plot.margin = ggplot2::margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                     axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm"),
          axis.title = element_text(face = "bold",size = 15),
          axis.text = element_text(size = 12))
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size = title.size,face = "bold"),
                                          axis.text.y = element_text(size = 12),
                                          axis.title.y = element_text(face = "bold",size = 15),)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- self.Description
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp,
                                                                         xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.6), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                          0.8))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = ggplot2::margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}




footer = function(){list(tags$footer(id = "myfooter",style="background-color:#033c73;position:fixed;height:120px;bottom:0px;width: 100%;right:0px;z-index:99999",
                                     fluidRow(
                                       
                                       div(img(src="TIMSource.png",style="width:120px;"),class="col-sm-3",style="text-align:right;"),
                                       div(class="col-sm-6",style="text-align:center;",
                                           a(href="https://www.pku-shenlab.cn",p("Shenlab, Peking University Cancer Hospital",style="color:white;font-weight:bolder;font-size:20px;text-align:center")),
                                           a(href="http://www.ncpsb.org/",p("Beijing Proteome Research Center, National Center for Protein Sciences (Beijing)",style="color:white;font-weight:bolder;font-size:20px;text-align:center")),
                                           a(href="http://www.lifelink.cn/sy",p("SIP LifeLink Oncology Research Institute",style="color:white;font-weight:bolder;font-size:20px;text-align:center"))
                                       ),
                                       div(class="col-sm-3",style="text-align:left;",
                                           html('<script type="text/javascript" src="//rf.revolvermaps.com/0/0/2.js?i=5z3og620gr1&amp;m=0&amp;s=120&amp;c=ff0000&amp;t=1" async="async"></script>')
                                       )
                                       
                                       
                                     )
                                     
                                     
                                     
))}

style_snvio <- function(){
  
    tags$style(
             HTML
             (
               "div#driver-popover-item{display:none;position:absolute;background:#fff;color:#2d2d2d;margin:0;padding:15px;border-radius:5px;min-width:250px;max-width:500px;box-shadow:0 1px 10px rgba(0,0,0,.4);z-index:1000000000;}",
               ".tooltip-inner{max-width: 500px;padding: 3px 8px;color: #ffffff;text-align: center;background-color: #000000;border-radius: 4px}",
               ".modal-footer{text-align:left}",
               "#second_menus1{background-color:#e7edf5;display:flex;flex-direction:row;justify-content: center;font-size:20px;margin-top:20px;margin-left:20px;margin-right:20px;border-radius:50px;border-block-color:#e7edf5;}",
               "#second_menus2{background-color:#e7edf5;display:flex;flex-direction:row;justify-content: center;font-size:20px;margin-top:20px;margin-left:20px;margin-right:20px;border-radius:50px;border-block-color:#e7edf5;}",
               "#second_menus3{background-color:#e7edf5;display:flex;flex-direction:row;justify-content: center;font-size:20px;margin-top:20px;margin-left:20px;margin-right:20px;border-radius:50px;border-block-color:#e7edf5;}",
               
               "#second_menus1 > li > a{border-radius: 40px;margin-top:2px;margin-bottom:2px;margin-right:50px;margin-left:50px;}",
               "#second_menus2 > li > a{border-radius: 40px;margin-top:2px;margin-bottom:2px;margin-right:50px;margin-left:50px;}",
               "#second_menus3 > li > a{border-radius: 40px;margin-top:2px;margin-bottom:2px;margin-right:50px;margin-left:50px;}",
               
               "#explore_ref_single{background-color:#033c73;font-size:20px;margin-top:20px;margin-left:0px;margin-right:0px;border-block-color:#033c73;}",
               "#explore_ref_single > li > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:white}",
               "#explore_ref_single > li.active > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:black}",
               
               "#explore_ref_pm{background-color:#033c73;font-size:20px;margin-top:20px;margin-left:0px;margin-right:0px;border-block-color:#033c73;}",
               "#explore_ref_pm > li > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:white}",
               "#explore_ref_pm > li.active > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:black}",
               
               "#explore_TCGA_single{background-color:#033c73;font-size:20px;margin-top:20px;margin-left:0px;margin-right:0px;border-block-color:#033c73;}",
               "#explore_TCGA_single > li > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:white}",
               "#explore_TCGA_single > li.active > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:black}",
               
               "#explore_TCGA_pm{background-color:#033c73;font-size:20px;margin-top:20px;margin-left:0px;margin-right:0px;border-block-color:#033c73;}",
               "#explore_TCGA_pm > li > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:white}",
               "#explore_TCGA_pm > li.active > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:black}",
               
               "#explore_CPTAC_single{background-color:#033c73;font-size:20px;margin-top:20px;margin-left:0px;margin-right:0px;border-block-color:#033c73;}",
               "#explore_CPTAC_single > li > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:white}",
               "#explore_CPTAC_single > li.active > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:black}",
               
               "#explore_CPTAC_pm{background-color:#033c73;font-size:20px;margin-top:20px;margin-left:0px;margin-right:0px;border-block-color:#033c73;}",
               "#explore_CPTAC_pm > li > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:white}",
               "#explore_CPTAC_pm > li.active > a{border-radius: 0px;margin-top:2px;margin-bottom:0px;color:black}",
               
               ".container-fluid{padding:0px}",
               '#fluid_id{padding:0px}',
               "#navbar_id{float:left;padding:0px;}",
               ".navbar{margin-bottom:0px;margin-left:-15px;margin-right:-15px}",
               ".row{margin-left:0px;margin-right:0px;}",
               ".skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper{background-color:white;}",
               ".content-wrapper{background-color:white;}",
               ".datatables{overflow-x:auto;}",
               
               "#Survival{text-align:center;overflow-x:auto;}",
               "#response{text-align:center;overflow-x:auto;}",
               "#TMB{text-align:center;overflow-x:auto;}",
               "#SPLOT{text-align:center;overflow-x:auto;}",
               "#immune_infiltration_datasets_all{text-align:center;overflow-x:auto;}",
               "#immune_signature_datasets_all{text-align:center;overflow-x:auto;}",
               "#immune_infiltration_datasets_one{text-align:center;overflow-x:auto;}",
               "#immune_signature_datasets_one{text-align:center;overflow-x:auto;}",
               
               "#Survival_pm{text-align:center;overflow-x:auto;}",
               "#response_pm{text-align:center;overflow-x:auto;}",
               "#TMB_pm{text-align:center;overflow-x:auto;}",
               "#SPLOT_pm{text-align:center;overflow-x:auto;}",
               "#immune_infiltration_datasets_all_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature_datasets_all_pm{text-align:center;overflow-x:auto;}",
               "#immune_infiltration_datasets_one_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature_datasets_one_pm{text-align:center;overflow-x:auto;}",
               
               "#immune_infiltration1{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2{text-align:center;overflow-x:auto;}",
               "#immune_signature1{text-align:center;overflow-x:auto;}",
               "#immune_signature2{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_tcga{text-align:center;overflow-x:auto;}",
               
               "#immune_infiltration1_pm{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature1_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature2_pm{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_tcga_pm{text-align:center;overflow-x:auto;}",
               
               "#immune_infiltration1_subtype{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_subtype{text-align:center;overflow-x:auto;}",
               "#immune_signature1_subtype{text-align:center;overflow-x:auto;}",
               "#immune_signature2_subtype{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_tcga_subtype{text-align:center;overflow-x:auto;}",
               
               "#immune_infiltration1_CPTAC_rna{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_CPTAC_rna{text-align:center;overflow-x:auto;}",
               "#immune_infiltration1_CPTAC_protein{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_CPTAC_protein{text-align:center;overflow-x:auto;}",
               "#immune_signature1_CPTAC_rna{text-align:center;overflow-x:auto;}",
               "#immune_signature2_CPTAC_rna{text-align:center;overflow-x:auto;}",
               "#immune_signature1_CPTAC_protein{text-align:center;overflow-x:auto;}",
               "#immune_signature2_CPTAC_protein{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_cptac_rna{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_cptac_protein{text-align:center;overflow-x:auto;}",
               
               "#immune_infiltration1_CPTAC_rna_pm{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_CPTAC_rna_pm{text-align:center;overflow-x:auto;}",
               "#immune_infiltration1_CPTAC_protein_pm{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_CPTAC_protein_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature1_CPTAC_rna_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature2_CPTAC_rna_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature1_CPTAC_protein_pm{text-align:center;overflow-x:auto;}",
               "#immune_signature2_CPTAC_protein_pm{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_cptac_rna_pm{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_cptac_protein_pm{text-align:center;overflow-x:auto;}",
               
               "#immune_infiltration1_CPTAC_rna_subtype{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_CPTAC_rna_subtype{text-align:center;overflow-x:auto;}",
               "#immune_infiltration1_CPTAC_protein_subtype{text-align:center;overflow-x:auto;}",
               "#immune_infiltration2_CPTAC_protein_subtype{text-align:center;overflow-x:auto;}",
               "#immune_signature1_CPTAC_rna_subtype{text-align:center;overflow-x:auto;}",
               "#immune_signature2_CPTAC_rna_subtype{text-align:center;overflow-x:auto;}",
               "#immune_signature1_CPTAC_protein_subtype{text-align:center;overflow-x:auto;}",
               "#immune_signature2_CPTAC_protein_subtype{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_cptac_rna_subtype{text-align:center;overflow-x:auto;}",
               "#GSEA_plot_cptac_protein_subtype{text-align:center;overflow-x:auto;}",
               
               "#volcano{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_pm{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_subtype{text-align:center;display:flex;align-items: center;justify-content: center;}",
               
               "#volcano_CPTAC_rna{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_CPTAC_protein{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_CPTAC_rna_pm{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_CPTAC_protein_pm{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_CPTAC_rna_subtype{text-align:center;display:flex;align-items: center;justify-content: center;}",
               "#volcano_CPTAC_protein_subtype{text-align:center;display:flex;align-items: center;justify-content: center;}",
               
               "#TCGA_survival{text-align:center;overflow-x:auto;}",
               "#TCGA_survival_pm{text-align:center;overflow-x:auto;}",
               "#TCGA_survival_subtype{text-align:center;overflow-x:auto;}",
               
               "#CPTAC_survival{text-align:center;overflow-x:auto;}",
               "#CPTAC_survival_pm{text-align:center;overflow-x:auto;}",
               "#CPTAC_survival_subtype{text-align:center;overflow-x:auto;}",

               "#sidebar_id1{color:white;background-color:#033c73;}",
               "#sidebar_id1_pm{color:white;background-color:#033c73;}",
               
               "#sidebar_id2{color:white;background-color:#033c73;}",
               "#sidebar_id2_pm{color:white;background-color:#033c73;}",
               "#sidebar_id2_subtype{color:white;background-color:#033c73;}",
               
               "#sidebar_id3{color:white;background-color:#033c73;}",
               "#sidebar_id3_pm{color:white;background-color:#033c73;}",
               "#sidebar_id3_subtype{color:white;background-color:#033c73;}",
               
               "//#warning{font-size:15px;font-weight:bolder;color:gray;text-align:center}",
               "//#warning_pm{font-size:15px;font-weight:bolder;color:gray;text-align:center}",
               "//#warning2{font-size:15px;font-weight:bolder;color:gray;text-align:center}",
               "//#warning2_pm{font-size:15px;font-weight:bolder;color:gray;text-align:center}",
               
               ".box-header{background-color:#033c73;}",
               ".box-title{color:white;}",
               ".box{border-top:3px solid #033c73;color}",
               ".btn-box-tool{color:white}",
               "#button {background-color: aliceblue;border-radius: 12px;color: #317eac;padding: 15px 32px;text-align: center;text-decoration: none;font-size: 16px; box-shadow: 0 8px 16px 0 rgba(0,0,0,0.2), 0 6px 20px 0 rgba(0,0,0,0.19)}",
               "#button2 {background-color: aliceblue;border-radius: 12px;color: #317eac;padding: 15px 32px;text-align: center;text-decoration: none;font-size: 16px; box-shadow: 0 8px 16px 0 rgba(0,0,0,0.2), 0 6px 20px 0 rgba(0,0,0,0.19)}"
               )
  )
}
