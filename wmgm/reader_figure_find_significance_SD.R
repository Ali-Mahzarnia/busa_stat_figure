library(ggplot2)

library("xlsx")
library("emmeans")
library("ggpubr")
library('openxlsx')

path_master = "AD_DECODE_data3.xlsx"
master = readxl::read_xlsx(path_master)
# master = t(na.omit(t(master)))



path_stats= "stats/"
stat_file_list_1=list.files(path_stats)
plain_index = grep("Tractometry_mrtrixfa", stat_file_list_1)
stat_file_list = stat_file_list_1[plain_index]
gmwm_index_file_list = stat_file_list_1[grep("Tractometry_greywhite_", stat_file_list_1)]


temp_data = read.csv(paste0(path_stats,stat_file_list[1]), sep = ",")
temp_data = temp_data [,2:dim(temp_data)[2]]
colnames = colnames(temp_data)
master = as.data.frame(master)
master = cbind(master, setNames( lapply(colnames, function(x) x=NA), colnames) )




for (i in 1:dim(master)[1]){
  
  index_masrter =  which(master$MRI_Exam[i]  == as.numeric(substr(stat_file_list,23,27)) )
  if (length(index_masrter) > 0){
    data = read.csv(paste0(path_stats,stat_file_list[index_masrter]), sep = ",")
    data = data[,2:dim(data)[2]]
    ####
    # dim(data)
    # c= c(seq(25:75))
    # data = data[c, ]
    ####
    ## taking into account grey and white
    index_gmwm = which(master$MRI_Exam[i] == as.numeric(substr(gmwm_index_file_list,24,28)) )
    data_gmwm = read.csv(paste0(path_stats,gmwm_index_file_list[index_gmwm]), sep = ",")
    data_gmwm = data_gmwm[,2:dim(data_gmwm)[2]]
    data_gmwm[data_gmwm==100] = 0
    data_gmwm[data_gmwm==101] = 1
    #commenting out or not the next line to not only look at white matter or the other way around
    data = data  * data_gmwm
    
    master[i,colnames] = sapply(data, sd)
    
    
  }
}  

geno = master$genotype
geno[master$genotype=="APOE34"] = "APOE44"
geno[master$genotype=="APOE23"] = "APOE33"
master$geno = geno

master = master[!is.na(master$genotype),]

master$age = as.numeric(master$age)

for (i in 48:(dim(master)[2]-1)) {
  temp = summary(lm(unlist(master[, i])~master$age* as.factor(master$geno ) ))
  temp = temp$coefficients
  if( !is.na(temp[4,4]) & temp[4,4] <= 0.05){
    print(i)
    
  plt = ggplot(master, aes(x=age, y= unlist(master[, i]) , color=geno, shape=geno)) +
    geom_point() +
    geom_smooth(method='lm', fullrange=TRUE, aes(fill = geno)) +
    theme_linedraw() +
    theme(text = element_text(size=20))+
    labs( title = paste0("P-Value =  ", temp[4,4]),  x = "Age",   y = paste0( "SD FA of " , colnames( master)[i]  ) ) 
  ggsave( paste0("figsSD/SD_",colnames( master)[ i], ".png") )
  }  else {print(temp[4,4])}
}






