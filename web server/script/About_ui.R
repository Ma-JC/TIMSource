About_ui = function(){
  tagList(
    br(),
    box(title = "Ref_ICI_datasets",width = 12,solidHeader = T,collapsible = T,collapsed = F,
        p("1. Samstein et al, Pan-cancer (Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Samstein, Robert M et al. 'Tumor mutational load predicts survival after immunotherapy across multiple cancer types.' Nature genetics vol. 51,2 (2019): 202-206. doi:10.1038/s41588-018-0312-8",
            href="https://www.nature.com/articles/s41588-018-0312-8"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),

        p("2. Van Allen et al, Melanoma (CTLA-4):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Van Allen, Eliezer M et al. 'Genomic correlates of response to CTLA-4 blockade in metastatic melanoma.' Science (New York, N.Y.) vol. 350,6257 (2015): 207-211. doi:10.1126/science.aad0095",
            href="https://www.science.org/doi/10.1126/science.aad0095"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),

        p("3. Miao et al, Microsatellite-stable solid tumors (Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Miao, Diana et al. 'Genomic correlates of response to immune checkpoint blockade in microsatellite-stable solid tumors.' Nature genetics vol. 50,9 (2018): 1271-1281. doi:10.1038/s41588-018-0200-2",
            href="https://www.nature.com/articles/s41588-018-0200-2"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("4. Janjigian et al, Esophagogastric Cancer (Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Janjigian, Yelena Y et al. 'Genetic Predictors of Response to Systemic Therapy in Esophagogastric Cancer.' Cancer discovery vol. 8,1 (2018): 49-58. doi:10.1158/2159-8290.CD-17-0787",
            href="https://cancerdiscovery.aacrjournals.org/content/8/1/49"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("5. Rizvi et al, NSCLC (Anti-PD-1/PD-L1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Rizvi, Hira et al. 'Molecular Determinants of Response to Anti-Programmed Cell Death (PD)-1 and Anti-Programmed Death-Ligand 1 (PD-L1) Blockade in Patients With Non-Small-Cell Lung Cancer Profiled With Targeted Next-Generation Sequencing.' Journal of clinical oncology : official journal of the American Society of Clinical Oncology vol. 36,7 (2018): 633-641. doi:10.1200/JCO.2017.75.3384",
            href="https://ascopubs.org/doi/10.1200/JCO.2017.75.3384"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("6. Pender, Pan-cancer (Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Pender, Alexandra et al. 'Genome and Transcriptome Biomarkers of Response to Immune Checkpoint Inhibitors in Advanced Solid Tumors.' Clinical cancer research : an official journal of the American Association for Cancer Research vol. 27,1 (2021): 202-212. doi:10.1158/1078-0432.CCR-20-1163",
            href="https://clincancerres.aacrjournals.org/content/27/1/202"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("7. Hellmann, et al. NSCLC (Anti-PD-1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Hellmann, Matthew D et al. 'Genomic Features of Response to Combination Immunotherapy in Patients with Advanced Non-Small-Cell Lung Cancer.' Cancer cell vol. 33,5 (2018): 843-852.e4. doi:10.1016/j.ccell.2018.03.018",
            href="https://www.sciencedirect.com/science/article/pii/S1535610818301235?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("8. Hugo, et al. Melanoma (Anti-PD1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Hugo, Willy et al. 'Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma.' Cell vol. 165,1 (2016): 35-44. doi:10.1016/j.cell.2016.02.065",
            href="https://www.sciencedirect.com/science/article/pii/S009286741630215X?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("9. Jiao, et al, Gastrointestinal cancer (Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Jiao, Xi et al. 'A genomic mutation signature predicts the clinical outcomes of immunotherapy and characterizes immunophenotypes in gastrointestinal cancer.' NPJ precision oncology vol. 5,1 36. 4 May. 2021, doi:10.1038/s41698-021-00172-5",
            href="https://www.nature.com/articles/s41698-021-00172-5"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("10. Snyder, et al. Melanoma (Anti-CTLA4):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Snyder, Alexandra et al. 'Genetic basis for clinical response to CTLA-4 blockade in melanoma.' The New England journal of medicine vol. 371,23 (2014): 2189-2199. doi:10.1056/NEJMoa1406498",
            href="https://www.nejm.org/doi/10.1056/NEJMoa1406498"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("11. Mariathasan, et al. Urothelial cancer (Anti-PDL1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Mariathasan, Sanjeev et al. 'TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells.' Nature vol. 554,7693 (2018): 544-548. doi:10.1038/nature25501",
            href="https://www.nature.com/articles/nature25501"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("12. Braun, et al. Clear cell renal cell carcinoma (Anti-PD1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Braun, David A et al. 'Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma.' Nature medicine vol. 26,6 (2020): 909-918. doi:10.1038/s41591-020-0839-y",
            href="https://www.nature.com/articles/s41591-020-0839-y"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("13. Motzer, et al, renal cell carcinoma (Anti-PD-L1+Axitinib):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Motzer, Robert J et al. 'Avelumab plus axitinib versus sunitinib in advanced renal cell carcinoma: biomarker analysis of the phase 3 JAVELIN Renal 101 trial.' Nature medicine vol. 26,11 (2020): 1733-1741. doi:10.1038/s41591-020-1044-8",
            href="https://www.nature.com/articles/s41591-020-1044-8"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("14. Miao et al, Clear cell renal cell carcinoma (Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Miao, Diana et al. 'Genomic correlates of response to immune checkpoint therapies in clear cell renal cell carcinoma.' Science (New York, N.Y.) vol. 359,6377 (2018): 801-806. doi:10.1126/science.aan5951",
            href="https://www.science.org/doi/10.1126/science.aan5951"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("15. Zhao et al. Glioblastoma(Anti-PD-1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Zhao, Junfei et al. 'Immune and genomic correlates of response to anti-PD-1 immunotherapy in glioblastoma.' Nature medicine vol. 25,3 (2019): 462-469. doi:10.1038/s41591-019-0349-y",
            href="https://www.nature.com/articles/s41591-019-0349-y"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("16. Rizvi et al.(2015) NSCLC(Anti-PD-1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Rizvi, Naiyer A et al. 'Cancer immunology. Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer.' Science (New York, N.Y.) vol. 348,6230 (2015): 124-8. doi:10.1126/science.aaa1348",
            href="https://www.science.org/doi/10.1126/science.aaa1348"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("17. Harding et al. Hepatocellular carcinoma(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Harding, James J et al. 'Prospective Genotyping of Hepatocellular Carcinoma: Clinical Implications of Next-Generation Sequencing for Matching Patients to Targeted and Immune Therapies.' Clinical cancer research : an official journal of the American Association for Cancer Research vol. 25,7 (2019): 2116-2126. doi:10.1158/1078-0432.CCR-18-2293",
            href="https://clincancerres.aacrjournals.org/content/25/7/2116"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("18. Anagnostou et al. Melanoma(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Anagnostou, Valsamo et al. 'Integrative Tumor and Immune Cell Multi-omic Analyses Predict Response to Immune Checkpoint Blockade in Melanoma.' Cell reports. Medicine vol. 1,8 100139. 17 Nov. 2020, doi:10.1016/j.xcrm.2020.100139",
            href="Integrative Tumor and Immune Cell Multi-omic Analyses Predict Response to Immune Checkpoint Blockade in Melanoma "),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("19. Anagnostou et al. NSCLC(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Anagnostou, Valsamo et al. 'Multimodal genomic features predict outcome of immune checkpoint blockade in non-small-cell lung cancer.' Nature cancer vol. 1,1 (2020): 99-111. doi:10.1038/s43018-019-0008-8",
            href="https://www.nature.com/articles/s43018-019-0008-8"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("20. Riaz et al. Melanoma(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Riaz, Nadeem et al. 'Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab.' Cell vol. 171,4 (2017): 934-949.e16. doi:10.1016/j.cell.2017.09.028",
            href="https://www.sciencedirect.com/science/article/pii/S0092867417311224?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("21. Liu et al. Melanoma(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Liu, David et al. 'Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma.' Nature medicine vol. 25,12 (2019): 1916-1927. doi:10.1038/s41591-019-0654-5",
            href="https://www.nature.com/articles/s41591-019-0654-5"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("22. Bai et al. Gastric cancer(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Bai, Yuezong et al. 'Efficacy and predictive biomarkers of immunotherapy in Epstein-Barr virus-associated gastric cancer.' Journal for immunotherapy of cancer vol. 10,3 (2022): e004080. doi:10.1136/jitc-2021-004080",
            href="https://jitc.bmj.com/content/10/3/e004080"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("23. Lu et al.  Neuroendocrine neoplasms(Anti-PD1/PDL1):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Lu, Ming et al. 'Efficacy, Safety, and Biomarkers of Toripalimab in Patients with Recurrent or Metastatic Neuroendocrine Neoplasms: A Multiple-Center Phase Ib Trial.' Clinical cancer research : an official journal of the American Association for Cancer Research vol. 26,10 (2020): 2337-2345. doi:10.1158/1078-0432.CCR-19-4000",
            href="https://aacrjournals.org/clincancerres/article/26/10/2337/82394/Efficacy-Safety-and-Biomarkers-of-Toripalimab-in"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("24. Wang et al. Gastrointestinal cancer(Mixed ICB):",style="font-weight: bolder;font-size: 16px;"),
        p(a("Wang, Zhenghang et al. 'Combination of AKT1 and CDH1 mutations predicts primary resistance to immunotherapy in dMMR/MSI-H gastrointestinal cancer.' Journal for immunotherapy of cancer vol. 10,6 (2022): e004703. doi:10.1136/jitc-2022-004703",
            href="https://jitc.bmj.com/content/10/6/e004703"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        ),
    
    
    box(title = "CPTAC_datanase",width = 12,solidHeader = T,collapsible = T,collapsed = F,
        p("1. Vasaikar et al. Colon cancer:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Vasaikar, Suhas et al. 'Proteogenomic Analysis of Human Colon Cancer Reveals New Therapeutic Opportunities.' Cell vol. 177,4 (2019): 1035-1049.e19. doi:10.1016/j.cell.2019.03.030",
            href="https://www.sciencedirect.com/science/article/pii/S0092867419302922?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("2. Krug et al. Breast cancer:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Krug, Karsten et al. 'Proteogenomic Landscape of Breast Cancer Tumorigenesis and Targeted Therapy.' Cell vol. 183,5 (2020): 1436-1456.e31. doi:10.1016/j.cell.2020.10.036",
            href="https://www.sciencedirect.com/science/article/pii/S0092867420314008?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("3. Clark et al. Clear cell renal cell carcinoma:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Clark, David J et al. 'Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma.' Cell vol. 179,4 (2019): 964-983.e31. doi:10.1016/j.cell.2019.10.007",
            href="https://www.sciencedirect.com/science/article/pii/S0092867419311237?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("4. Huang et al. HPV-negative head and neck squamous cell carcinomas:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Huang, Chen et al. 'Proteogenomic insights into the biology and treatment of HPV-negative head and neck squamous cell carcinoma.' Cancer cell vol. 39,3 (2021): 361-379.e16. doi:10.1016/j.ccell.2020.12.007",
            href="https://www.sciencedirect.com/science/article/pii/S1535610820306553?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("5. Satpathy et al. Lung squamous cell carcinoma:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Satpathy, Shankha et al. 'A proteogenomic portrait of lung squamous cell carcinoma.' Cell vol. 184,16 (2021): 4348-4371.e40. doi:10.1016/j.cell.2021.07.016",
            href="https://www.sciencedirect.com/science/article/pii/S0092867421008576?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("6. Gillette et al. Lung adenocarcinomas:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Gillette, Michael A et al. 'Proteogenomic Characterization Reveals Therapeutic Vulnerabilities in Lung Adenocarcinoma.' Cell vol. 182,1 (2020): 200-225.e35. doi:10.1016/j.cell.2020.06.013",
            href="https://www.sciencedirect.com/science/article/pii/S0092867420307443?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("7. Hu et al & McDermott et al. High-grade serous ovarian carcinomas:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Hu, Yingwei et al. 'Integrated Proteomic and Glycoproteomic Characterization of Human High-Grade Serous Ovarian Carcinoma.' Cell reports vol. 33,3 (2020): 108276. doi:10.1016/j.celrep.2020.108276",
            href="https://www.sciencedirect.com/science/article/pii/S2211124720312651?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        p(a("McDermott, Jason E et al. 'Proteogenomic Characterization of Ovarian HGSC Implicates Mitotic Kinases, Replication Stress in Observed Chromosomal Instability.' Cell reports. Medicine vol. 1,1 (2020): 100004. doi:10.1016/j.xcrm.2020.100004",
            href="https://www.sciencedirect.com/science/article/pii/S2666379120300045?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("8. Dou et al. Uterine Corpus Endometrial Carcinoma:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Dou, Yongchao et al. 'Proteogenomic Characterization of Endometrial Carcinoma.' Cell vol. 180,4 (2020): 729-748.e26. doi:10.1016/j.cell.2020.01.026",
            href="https://www.sciencedirect.com/science/article/pii/S0092867420301070?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("9. Wang et al.Glioblastoma:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Wang, Liang-Bo et al. 'Proteogenomic and metabolomic characterization of human glioblastoma.' Cancer cell vol. 39,4 (2021): 509-528.e20. doi:10.1016/j.ccell.2021.01.006",
            href="https://www.sciencedirect.com/science/article/pii/S1535610821000507?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("10. Cao et al.pancreatic ductal adenocarcinoma:",style="font-weight: bolder;font-size: 16px;"),
        p(a("Cao, Liwei et al. 'Proteogenomic characterization of pancreatic ductal adenocarcinoma.' Cell vol. 184,19 (2021): 5031-5052.e26. doi:10.1016/j.cell.2021.08.023",
            href="https://www.sciencedirect.com/science/article/pii/S0092867421009971?via%3Dihub"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        
        ),
    box(title = "TCGA_database",width = 12,solidHeader = T,collapsible = T,collapsed = F,
        
        p("PanCanAtlas | RNA:",style="font-weight: bolder;font-size: 16px;"),
        p(a("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",
            href="http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("PanCanAtlas | Mutations:",style="font-weight: bolder;font-size: 16px;"),
        p(a("mc3.v0.2.8.PUBLIC.maf.gz",
            href="http://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        p("PanCanAtlas | Clinical Data Resource(CDR):",style="font-weight: bolder;font-size: 16px;"),
        p(a("TCGA-CDR-SupplementalTableS1.xlsx",
            href="https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"),
          style="text-indent:1cm;font-size:15px;font-weight:bold;"),
        
        
        ),
    div(style="height:200px;display: inline-block;width: 20px;")

  )
}