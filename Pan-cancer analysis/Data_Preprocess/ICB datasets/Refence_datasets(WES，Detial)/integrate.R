dataset3_1.1 = read.table("./dataset3/result_data/dataset3_Bladder_Cancer.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset3_1.2 = read.table("./dataset3/result_data/dataset3_Melanoma.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset3_1.3 = read.table("./dataset3/result_data/dataset3_Non-Small_Cell_Lung_Cancer.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)


dataset2_1 = read.table("./dataset2/result_data/dataset2.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset3_1 = read.table("./dataset3/result_data/dataset3.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset6_1 = read.table("./dataset6/result_data/dataset6.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset7_1 = read.table("./dataset7/result_data/dataset7.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset8_1 = read.table("./dataset8/result_data/dataset8.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset9_1 = read.table("./dataset9/result_data/dataset9.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset11_1 = read.table("./dataset10/result_data/dataset10.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset12_1 = read.table("./dataset12/result_data/dataset12.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset13_1 = read.table("./dataset13/result_data/dataset13.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset14_1 = read.table("./dataset14/result_data/dataset14.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset15_1 = read.table("./dataset15/result_data/dataset15.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset17_1 = read.table("./dataset16/result_data/dataset16.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset18_1 = read.table("./dataset18/result_data/dataset18.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset19_1 = read.table("./dataset19/result_data/dataset19.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset20_1 = read.table("./dataset20/result_data/dataset20.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
dataset22_1 = read.table("./dataset21/result_data/dataset21.txt",sep = "\t",quote = "",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

data = list(

  "dataset2" = dataset2_1,
  "dataset3.1" = dataset3_1.1,
  "dataset3.2" = dataset3_1.2,
  "dataset3.3" = dataset3_1.3,
  "dataset3" = dataset3_1,
  "dataset6" = dataset6_1,
  "dataset7" = dataset7_1,
  "dataset8" = dataset8_1,
  "dataset9" = dataset9_1,
  "dataset10" = dataset11_1,
  "dataset12" = dataset12_1,
  "dataset13" = dataset13_1,
  "dataset14" = dataset14_1,
  "dataset15" = dataset15_1,
  "dataset16" = dataset17_1,
  "dataset18" = dataset18_1,
  "dataset19" = dataset19_1,
  "dataset20" = dataset20_1,
  "dataset21" = dataset22_1,
)

colnames(data$dataset7)[7] = "origin_therapy"
data$dataset7$origin_therapy = "Not Available"


data$dataset16$origin_therapy = "Not Available"



sapply(data,nrow)

saveRDS(data,"Refence_datasets_COX.rds")
