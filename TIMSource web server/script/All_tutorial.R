#### ref + single####
guide <- Cicerone$
  new()$
  step(
    el = "second_menus1",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, you have chosen immunotherapy-related datasets. For further details of the datasets please refer to the 'About' section."
  )$
  step(
    el = "sidebar_id1",
    title = "Help",
    description = "This section allows the users to select the dataset or cancer type, the mutation gene or pathway for patient stratification, and the mutation type."
  )$
  step(
    el = "dataset_tutorial",
    title = "Help",
    description = tooltip_text["dataset",]
  )$
  step(
    "gene_symbol_tutorial",
    "Help",
    tooltip_text["gene_symbol",]
  )$
  step(
    "Mut_type_ref_single_tutorial",
    "Help",
    tooltip_text["Variant_Classification",]
  )$
  step(
    "Wild_type_ref_single_tutorial",
    "Help",
    tooltip_text["Wildtype groups",]
  )$
  step(
    "Therapy_type_ref_single_tutorial",
    "Help",
    tooltip_text["Immune_checkpoint_blockade",]
  )$
  step(
    "sec",
    "Help",
    "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "ref_single_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id1",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )

#### tcga + single ####
guide_tcga_single <- Cicerone$
  new()$
  step(
    el = "second_menus1",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen TCGA database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id2",
    title = "Help",
    description = "This section allows the users to select the dataset or cancer type, the mutation gene or pathway for patient stratification, and the mutation type."
  )$
  step(
    el = "Cancer_type_tutorial",
    title = "Help",
    description = tooltip_text["cancer type tcga",]
  )$
  step(
    el = "gene_symbol_TCGA_tutorial",
    title = "Help",
    description = tooltip_text["gene_symbol_tcga",]
  )$
  step(
    el = "Mut_type_tcga_single_tutorial",
    title = "Help",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    el = "Wild_type_tcga_single_tutorial",
    title = "Help",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    el = "sec2",
    title = "Help",
    description = "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "tcga_single_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id2",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )


#### cptac + single ####
guide_cptac_single <- Cicerone$
  new()$
  step(
    el = "second_menus1",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen CPTAC database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id3",
    title = "Help",
    description = "This section allows the users to select the dataset or cancer type, the mutation gene or pathway for patient stratification, and the mutation type."
  )$
  step(
    el = "Cancer_type2_tutorial",
    title = "Help",
    description = tooltip_text["cancer type cptac",]
  )$
  step(
    el = "gene_symbol_CPTAC_tutorial",
    title = "Help",
    description = tooltip_text["gene_symbol_tcga",]
  )$
  step(
    el = "Mut_type_cptac_single_tutorial",
    title = "Help",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    el = "Wild_type_cptac_single_tutorial",
    title = "Help",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    el = "sec3",
    title = "Help",
    description = "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "cptac_single_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id3",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )


#### ref + pathway ####
guide_pm <- Cicerone$
  new()$
  step(
    el = "second_menus2",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen immunotherapy-related datasets. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id1_pm",
    title = "Help",
    description = "This section allows the users to select the dataset or cancer type, the mutation gene or pathway for patient stratification, and the mutation type."
  )$
  step(
    el = "dataset_pm_tutorial",
    title = "Help",
    description = tooltip_text["dataset",]
  )$
  step(
    el = "Keyword_tutorial",
    title = "Help",
    description = tooltip_text["keyword",]
  )$
  step(
    el = "SUB_WORD_tutorial",
    title = "Help",
    description = "When you enter keywords, you need to click this button. At this point, the backend will search for relevant functions or pathways in MSigDB and display them in the selection box below."
  )$
  step(
    "gene_symbol_pm_tutorial",
    "Help",
    tooltip_text["pathway",]
  )$
  step(
    "Mut_type_ref_pm_tutorial",
    "Help",
    tooltip_text["Variant_Classification",]
  )$
  step(
    "Wild_type_ref_pm_tutorial",
    "Help",
    tooltip_text["Wildtype groups",]
  )$
  step(
    "Therapy_type_ref_pm_tutorial",
    "Help",
    tooltip_text["Immune_checkpoint_blockade",]
  )$
  step(
    "sec_pm",
    "Help",
    "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "ref_pm_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id1_pm",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )

#### TCGA + pathway ####
guide_tcga_pm <- Cicerone$
  new()$
  step(
    el = "second_menus2",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen TCGA database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id2_pm",
    title = "Help",
    description = "This section allows the users to select the dataset or cancer type, the mutation gene or pathway for patient stratification, and the mutation type."
  )$
  step(
    el = "Cancer_type_pm_tutorial",
    title = "Help",
    description = tooltip_text["cancer type tcga",]
  )$
  step(
    el = "Keyword_tcga_tutorial",
    title = "Help",
    description = tooltip_text["keyword",]
  )$
  step(
    el = "SUB_WORD_tcga_tutorial",
    title = "Help",
    description = "When you enter keywords, you need to click this button. At this point, the backend will search for relevant functions or pathways in MSigDB and display them in the selection box below."
  )$
  step(
    "gene_symbol_TCGA_pm_tutorial",
    "Help",
    description = tooltip_text["pathway",]
  )$
  step(
    "Mut_type_tcga_pm_tutorial",
    "Help",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_tcga_pm_tutorial",
    "Help",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "sec2_pm",
    "Help",
    description = "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "tcga_pm_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id2_pm",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )

#### cptac + pathway ####
guide_cptac_pm <- Cicerone$
  new()$
  step(
    el = "second_menus2",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen CPTAC database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id3_pm",
    title = "Help",
    description = "This section allows the users to select the dataset or cancer type, the mutation gene or pathway for patient stratification, and the mutation type."
  )$
  step(
    el = "Cancer_type3_pm_tutorial",
    title = "Help",
    description = tooltip_text["cancer type cptac",]
  )$
  step(
    el = "Keyword_cptac_tutorial",
    title = "Help",
    description = tooltip_text["keyword",]
  )$
  step(
    el = "SUB_WORD_cptac_tutorial",
    title = "Help",
    description = "When you enter keywords, you need to click this button. At this point, the backend will search for relevant functions or pathways in MSigDB and display them in the selection box below."
  )$
  step(
    "gene_symbol_CPTAC_pm_tutorial",
    "Help",
    description = tooltip_text["pathway",]
  )$
  step(
    "Mut_type_cptac_pm_tutorial",
    "Help",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_cptac_pm_tutorial",
    "Help",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "sec3_pm",
    "Help",
    description = "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "cptac_pm_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id3_pm",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )

#### tcga + subtype ####
guide_tcga_subtype <- Cicerone$
  new()$
  step(
    el = "second_menus3",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen TCGA database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id2_subtype",
    title = "Help",
    description = "This parameter bar is used to confirm the tumor type we want to analyze in detail, the gene mutation background, the definition of the mutation group and wildtype group, and the gene used to stratify patients based on expression levels."
  )$
  step(
    el = "Cancer_type_subtype_tutorial",
    title = "Help",
    description = tooltip_text["cancer type tcga",]
  )$
  step(
    el = "Subtype_tcga_tutorial",
    title = "Help",
    description = tooltip_text["WORM",]
  )$
  # step(
  #   el = "MT_OR_WT_tcga_tutorial",
  #   title = "Help",
  #   description = "This is where you enter the text you want to print."
  # )$
  step(
    "Mut_type_tcga_subtype_tutorial",
    "Help",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_tcga_subtype_tutorial",
    "Help",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "gene_symbol_TCGA_subtype_tutorial",
    "Help",
    description = tooltip_text["Gene_Symbol",]
  )$
  step(
    "sec2_subtype",
    "Help",
    description = "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "tcga_subtype_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id2_subtype",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )

#### cptac + subtype ####
guide_cptac_subtype <- Cicerone$
  new()$
  step(
    el = "second_menus3",
    title = "Help",
    description = "We provide three major data sources to explore the impact of the gene mutation or altered pathway on ICB outcomes and tumor immune characteristics in individual cohort or cancer type. Here, we have chosen CPTAC database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id3_subtype",
    title = "Help",
    description = "This parameter bar is used to confirm the tumor type we want to analyze in detail, the gene mutation background, the definition of the mutation group and wildtype group, and the gene used to stratify patients based on expression levels."
  )$
  step(
    el = "Cancer_type2_subtype_tutorial",
    title = "Help",
    description = tooltip_text["cancer type cptac",]
  )$
  step(
    el = "Subtype_cptac_tutorial",
    title = "Help",
    description = tooltip_text["WORM",]
  )$
  # step(
  #   el = "MT_OR_WT_cptac_tutorial",
  #   title = "Help",
  #   description = "This is where you enter the text you want to print."
  # )$
  step(
    "Mut_type_cptac_subtype_tutorial",
    "Help",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_cptac_subtype_tutorial",
    "Help",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "gene_symbol_CPTAC_subtype_tutorial",
    "Help",
    description = tooltip_text["Gene_Symbol",]
  )$
  # step(
  #   "gene_symbol_CPTAC_subtype2_tutorial",
  #   "Help",
  #   "Help to the server for printing"
  # )$
  step(
    el = "sec3_subtype",
    title = "Help",
    description = "When you have confirmed that all the selected parameters, you can click the 'Submit' button to initiate the calculation process. Once submission, the backend will perform the calculation, and the results will be displayed on the right side. Please note that this button cannot be clicked again until the calculation has finished."
  )$
  step(
    el = "cptac_subtype_analysis",
    title = "Help",
    description = "This tab displays the analytical functions that can be performed under the current data source. In ICB, TCGA and CPTAC dataset, users can obtain the mutation landscape of one gene or pathway, and perform the K-M survival analysis (TCGA and CPTAC can be considered as nonimmunotherapy treatment), response comparison (only in ICB dataset) according to the mutation status.Other General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis, which can be performed in TCGA and CPTAC dataset, and some of the ICB cohorts (if the transcriptional data was available in the published clinical study)."
  )$
  step(
    el = "main_id3_subtype",
    title = "Help",
    description = "Here, you will find all the results of your analysis. In order to improve website response time, some functions are not performed by default, and can be initiated by clicking the '+' button. Let's start exploring now!"
  )

