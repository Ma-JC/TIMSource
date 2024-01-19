#### ref + single####
guide <- Cicerone$
  new()$
  step(
    el = "second_menus1",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen immunotherapy-related datasets. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id1",
    title = "Text Input",
    description = "This parameter bar is used to confirm the dataset or cancer type we want to analyze in detail, the mutation gene or pathway we want to use to stratify patients, and how to define the mutation group and wildtype group."
  )$
  step(
    el = "dataset_tutorial",
    title = "Text Input",
    description = tooltip_text["dataset",]
  )$
  step(
    "gene_symbol_tutorial",
    "Send the Text",
    tooltip_text["gene_symbol",]
  )$
  step(
    "Mut_type_ref_single_tutorial",
    "Send the Text",
    tooltip_text["Variant_Classification",]
  )$
  step(
    "Wild_type_ref_single_tutorial",
    "Send the Text",
    tooltip_text["Wildtype groups",]
  )$
  step(
    "sec",
    "Send the Text",
    "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "ref_single_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id1",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )

#### tcga + single ####
guide_tcga_single <- Cicerone$
  new()$
  step(
    el = "second_menus1",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen TCGA database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id2",
    title = "Text Input",
    description = "This parameter bar is used to confirm the dataset or cancer type we want to analyze in detail, the mutation gene or pathway we want to use to stratify patients, and how to define the mutation group and wildtype group."
  )$
  step(
    el = "Cancer_type_tutorial",
    title = "Text Input",
    description = tooltip_text["cancer type tcga",]
  )$
  step(
    el = "gene_symbol_TCGA_tutorial",
    title = "Text Input",
    description = tooltip_text["gene_symbol_tcga",]
  )$
  step(
    el = "Mut_type_tcga_single_tutorial",
    title = "Text Input",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    el = "Wild_type_tcga_single_tutorial",
    title = "Text Input",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    el = "sec2",
    title = "Text Input",
    description = "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "tcga_single_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id2",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )


#### cptac + single ####
guide_cptac_single <- Cicerone$
  new()$
  step(
    el = "second_menus1",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen CPTAC database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id3",
    title = "Text Input",
    description = "This parameter bar is used to confirm the dataset or cancer type we want to analyze in detail, the mutation gene or pathway we want to use to stratify patients, and how to define the mutation group and wildtype group."
  )$
  step(
    el = "Cancer_type2_tutorial",
    title = "Text Input",
    description = tooltip_text["cancer type cptac",]
  )$
  step(
    el = "gene_symbol_CPTAC_tutorial",
    title = "Text Input",
    description = tooltip_text["gene_symbol_tcga",]
  )$
  step(
    el = "Mut_type_cptac_single_tutorial",
    title = "Text Input",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    el = "Wild_type_cptac_single_tutorial",
    title = "Text Input",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    el = "sec3",
    title = "Text Input",
    description = "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "cptac_single_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id3",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )


#### ref + pathway ####
guide_pm <- Cicerone$
  new()$
  step(
    el = "second_menus2",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen immunotherapy-related datasets. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id1_pm",
    title = "Text Input",
    description = "This parameter bar is used to confirm the dataset or cancer type we want to analyze in detail, the mutation gene or pathway we want to use to stratify patients, and how to define the mutation group and wildtype group."
  )$
  step(
    el = "dataset_pm_tutorial",
    title = "Text Input",
    description = tooltip_text["dataset",]
  )$
  step(
    el = "Keyword_tutorial",
    title = "Text Input",
    description = tooltip_text["keyword",]
  )$
  step(
    el = "SUB_WORD_tutorial",
    title = "Text Input",
    description = "When you enter keywords, you need to click this button. At this point, the backend will search for relevant functions or pathways in MSigDB and display them in the selection box below."
  )$
  step(
    "gene_symbol_pm_tutorial",
    "Send the Text",
    tooltip_text["pathway",]
  )$
  step(
    "Mut_type_ref_pm_tutorial",
    "Send the Text",
    tooltip_text["Variant_Classification",]
  )$
  step(
    "Wild_type_ref_pm_tutorial",
    "Send the Text",
    tooltip_text["Wildtype groups",]
  )$
  step(
    "sec_pm",
    "Send the Text",
    "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "ref_pm_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id1_pm",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )

#### TCGA + pathway ####
guide_tcga_pm <- Cicerone$
  new()$
  step(
    el = "second_menus2",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen TCGA database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id2_pm",
    title = "Text Input",
    description = "This parameter bar is used to confirm the dataset or cancer type we want to analyze in detail, the mutation gene or pathway we want to use to stratify patients, and how to define the mutation group and wildtype group."
  )$
  step(
    el = "Cancer_type_pm_tutorial",
    title = "Text Input",
    description = tooltip_text["cancer type tcga",]
  )$
  step(
    el = "Keyword_tcga_tutorial",
    title = "Text Input",
    description = tooltip_text["keyword",]
  )$
  step(
    el = "SUB_WORD_tcga_tutorial",
    title = "Text Input",
    description = "When you enter keywords, you need to click this button. At this point, the backend will search for relevant functions or pathways in MSigDB and display them in the selection box below."
  )$
  step(
    "gene_symbol_TCGA_pm_tutorial",
    "Send the Text",
    description = tooltip_text["pathway",]
  )$
  step(
    "Mut_type_tcga_pm_tutorial",
    "Send the Text",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_tcga_pm_tutorial",
    "Send the Text",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "sec2_pm",
    "Send the Text",
    description = "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "tcga_pm_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id2_pm",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )

#### cptac + pathway ####
guide_cptac_pm <- Cicerone$
  new()$
  step(
    el = "second_menus2",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen CPTAC database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id3_pm",
    title = "Text Input",
    description = "This parameter bar is used to confirm the dataset or cancer type we want to analyze in detail, the mutation gene or pathway we want to use to stratify patients, and how to define the mutation group and wildtype group."
  )$
  step(
    el = "Cancer_type3_pm_tutorial",
    title = "Text Input",
    description = tooltip_text["cancer type cptac",]
  )$
  step(
    el = "Keyword_cptac_tutorial",
    title = "Text Input",
    description = tooltip_text["keyword",]
  )$
  step(
    el = "SUB_WORD_cptac_tutorial",
    title = "Text Input",
    description = "When you enter keywords, you need to click this button. At this point, the backend will search for relevant functions or pathways in MSigDB and display them in the selection box below."
  )$
  step(
    "gene_symbol_CPTAC_pm_tutorial",
    "Send the Text",
    description = tooltip_text["pathway",]
  )$
  step(
    "Mut_type_cptac_pm_tutorial",
    "Send the Text",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_cptac_pm_tutorial",
    "Send the Text",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "sec3_pm",
    "Send the Text",
    description = "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "cptac_pm_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id3_pm",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )

#### tcga + subtype ####
guide_tcga_subtype <- Cicerone$
  new()$
  step(
    el = "second_menus3",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen TCGA database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id2_subtype",
    title = "Text Input",
    description = "This parameter bar is used to confirm the tumor type we want to analyze in detail, the gene mutation background, the definition of the mutation group and wildtype group, and the gene used to stratify patients based on expression levels."
  )$
  step(
    el = "Cancer_type_subtype_tutorial",
    title = "Text Input",
    description = tooltip_text["cancer type tcga",]
  )$
  step(
    el = "Subtype_tcga_tutorial",
    title = "Text Input",
    description = tooltip_text["WORM",]
  )$
  # step(
  #   el = "MT_OR_WT_tcga_tutorial",
  #   title = "Text Input",
  #   description = "This is where you enter the text you want to print."
  # )$
  step(
    "Mut_type_tcga_subtype_tutorial",
    "Send the Text",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_tcga_subtype_tutorial",
    "Send the Text",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "gene_symbol_TCGA_subtype_tutorial",
    "Send the Text",
    description = tooltip_text["Gene_Symbol",]
  )$
  step(
    "sec2_subtype",
    "Send the Text",
    description = "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "tcga_subtype_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id2_subtype",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )

#### cptac + subtype ####
guide_cptac_subtype <- Cicerone$
  new()$
  step(
    el = "second_menus3",
    title = "Text Input",
    description = "We provide three major data sources to analyze the relationship between gene mutations or pathway mutations and immunotherapy and the immune microenvironment. Here, we have chosen CPTAC database. For more information on the data, please see the 'About' section."
  )$
  step(
    el = "sidebar_id3_subtype",
    title = "Text Input",
    description = "This parameter bar is used to confirm the tumor type we want to analyze in detail, the gene mutation background, the definition of the mutation group and wildtype group, and the gene used to stratify patients based on expression levels."
  )$
  step(
    el = "Cancer_type2_subtype_tutorial",
    title = "Text Input",
    description = tooltip_text["cancer type cptac",]
  )$
  step(
    el = "Subtype_cptac_tutorial",
    title = "Text Input",
    description = tooltip_text["WORM",]
  )$
  # step(
  #   el = "MT_OR_WT_cptac_tutorial",
  #   title = "Text Input",
  #   description = "This is where you enter the text you want to print."
  # )$
  step(
    "Mut_type_cptac_subtype_tutorial",
    "Send the Text",
    description = tooltip_text["Variant_Classification_tcga",]
  )$
  step(
    "Wild_type_cptac_subtype_tutorial",
    "Send the Text",
    description = tooltip_text["Wildtype groups",]
  )$
  step(
    "gene_symbol_CPTAC_subtype_tutorial",
    "Send the Text",
    description = tooltip_text["Gene_Symbol",]
  )$
  # step(
  #   "gene_symbol_CPTAC_subtype2_tutorial",
  #   "Send the Text",
  #   "Send the text to the server for printing"
  # )$
  step(
    el = "sec3_subtype",
    title = "Text Input",
    description = "When you have confirmed that all the selected parameters are correct, you can click this button to submit. After the submission, the backend will perform the calculation based on the content involved in the current tab position on the right. At this point, this button cannot be clicked again until the calculation is finished."
  )$
  step(
    el = "cptac_subtype_analysis",
    title = "Text Input",
    description = "This tab bar displays the analysis functions that can be performed under the currently selected data source, which can be divided into general functions and functions unique to the data source. General functions mainly include immune infiltration, immune signature, differential expression analysis, and GSEA analysis. It is worth noting that only some of the Ref_ICI datasets with RNA-seq data can perform general functions. Although survival analysis is a function that can be performed by all three data sources, Ref_ICI datasets provide the patient survival situation under immunotherapy, which is different from TCGA and CPTAC. In addition, the unique functions of Ref_ICI datasets include analyzing the response to immune checkpoint inhibitors. TCGA and CPTAC data have more standardized and unified mutation data, so they can provide more detailed analysis at the mutation level (Ref_ICI datasets come from different studies, so only simple mutation types and mutation details for each patient are provided). CPTAC, because of its unique proteomics data, can also perform general functions at the protein level."
  )$
  step(
    el = "main_id3_subtype",
    title = "Text Input",
    description = "All analysis results will be displayed here. In order to improve the response speed of the web page, the results of some functions are hidden by default, and only by clicking the '+' button, the box will expand and the relevant analysis will be launched. Let's start exploring now!"
  )

