library(ggplot2)

path_master = "AD_DECODE_data3.xlsx"
master = readxl::read_xlsx(path_master)
# master = t(na.omit(t(master)))



path_stats= "stats/"
stat_file_list_1=list.files(path_stats)
plain_index = grep("Tractometry_mrtrixfa", stat_file_list_1)
stat_file_list = stat_file_list_1[plain_index]
gmwm_index_file_list = stat_file_list_1[grep("Tractometry_greywhite_", stat_file_list_1)]


temp_data = read.csv(paste0(path_stats,stat_file_list[1]), sep = ",")
colnames = colnames(temp_data)
colnames_index_left =grep( "left", colnames)
colnames_left = colnames[colnames_index_left]
colnames_no_leftright = gsub ("_left" , "", colnames_left)


master = as.data.frame(master)
master = cbind(master, setNames( lapply(colnames_no_leftright, function(x) x=0), colnames_no_leftright) )


for (i in 1:dim(master)[1]){
  
  index_masrter =  which(master$MRI_Exam[i]  == as.numeric(substr(stat_file_list,23,27))  )
  if (length(index_masrter) > 0){
  data = read.csv(paste0(path_stats,stat_file_list[index_masrter]), sep = ",")
  data = data[,2:dim(data)[2]]
  ## taking into account grey and white
  index_gmwm = which(master$MRI_Exam[i] == as.numeric(substr(gmwm_index_file_list,24,28)) )
  data_gmwm = read.csv(paste0(path_stats,gmwm_index_file_list[index_gmwm]), sep = ",")
  data_gmwm = data_gmwm[,2:dim(data_gmwm)[2]]
  data_gmwm[data_gmwm==100] = 0
  data_gmwm[data_gmwm==101] = 1
  #commenting out or not the next line to not only look at white matter or the other way around
  data = data  * data_gmwm
  
  
  for (j in colnames_left) {
    
    L = data[,j]
    R = data[, gsub("left", "right",j )]
    master[i,gsub("_left", "" ,j )] = mean( abs(L - R) / (L+R ) )
    # j=colnames_left[1]
    # master[i,gsub("_left", "" ,j )] =mean(data[,j] - data[, gsub("left", "right",j )] / (data[,j] + data[, gsub("left", "right",j )]))
    # master[i,gsub("_left", "" ,j )] =abs(mean(data[,j]) - mean(data[, gsub("left", "right",j )])) / mean(data[,j] + data[, gsub("left", "right",j )])
    # master[i,gsub("_left", "" ,j )] =mean(abs(data[,j] - data[, gsub("left", "right",j )]) / (data[,j] + data[, gsub("left", "right",j )]))    
    # master[i,gsub("_left", "" ,j )] =abs(mean(data[,j] - data[, gsub("left", "right",j )] / (data[,j] + data[, gsub("left", "right",j )])))
    
    }
  ####
  # dim(data)
  # c= c(seq(25:75))
  # data = data[c, ]
  ####
  # master[i,colnames] = sapply(data, mean)

  
  }
}  

geno = master$genotype
geno[master$genotype=="APOE34"] = "APOE44"
geno[master$genotype=="APOE23"] = "APOE33"
master$geno = geno

master = master[!is.na(master$genotype),]


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
    labs( title = paste0("P-Value = ", temp[4,4]),  x = "Age",   y = paste0( "AI of FA for " , colnames( master)[i]  ) ) 
  ggsave( paste0("figsAI/AI_",colnames( master)[ i], ".png") )
  }  else {print(temp[4,4])}
}






