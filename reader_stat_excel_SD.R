library(ggplot2)

library("xlsx")
library("emmeans")
library("ggpubr")
library('openxlsx')

path_master = "AD_DECODE_data3.xlsx"
master = readxl::read_xlsx(path_master)
# master = t(na.omit(t(master)))



path_stats= "stats/"
stat_file_list=list.files(path_stats)
plain_index = grep("Tractometry_mrtrixfa", stat_file_list)
stat_file_list = stat_file_list[plain_index]


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
    master[i,colnames] = sapply(data, sd)
    
    
  }
}  

geno = master$genotype
geno[master$genotype=="APOE34"] = "APOE44"
geno[master$genotype=="APOE23"] = "APOE33"
master$geno = geno

master = master[!is.na(master$genotype),]


##saving master for three age group#
masterorig = master

colnames(master)[dim(master)[2]] = "Genotype"
master$Genotype = as.factor(master$Genotype)
master$age= as.numeric(master$age)



for (i in 48:(dim(master)[2]-1)) {

  lm = lm(unlist(master[, i])~  age *  Genotype  , data=master)

  write.xlsx2(as.data.frame(  anova(lm)), paste0("sd.xlsx"), sheet =colnames(master)[i],append = TRUE)
  eta=effectsize::cohens_f(lm, alternative='two.sided')
  write.xlsx2(eta, paste0("sd.xlsx"), sheet =paste0(colnames(master)[i],"_cohen_"),append = TRUE)
  # res<-emmeans(lm, list(pairwise ~ factor(geno)), adjust="tukey")
  posthoc2 <- na.omit((emmeans(lm, list(pairwise ~ Genotype), adjust="tukey")$`pairwise differences of Genotype`)) 
  write.xlsx2(posthoc2, paste0("sd.xlsx"), sheet =paste0(colnames(master)[i],"_emmeans_"),append = TRUE)
  
  # contrast(res[[1]], "eff", by = "geno")
}



# young_thresh = quantile( master$age , 0.33)
# midage_thresh = quantile( master$age , 0.66)
# 
# 
# 
# 
# ######## first find out which bundle is significant at least in one of the three age groups 
# master = masterorig[masterorig$age <= young_thresh,]
# 
# Young_sig = matrix (NA, length(49:(dim(master)[2]-1)))
# for (i in 49:(dim(master)[2]-1)) {
#   temp = summary(lm(unlist(master[, i])~master$age + master$age: as.factor(master$geno ) ))
#   temp = temp$coefficients
#   if (temp[3,4]<0.05){
#     Young_sig[(i-48)] = "Sig"
#     # print(temp[3,4])
#   }
#   else { Young_sig[(i-48)] = "Non_sig" }
# }
# 
# 
# 
# master = masterorig[ (masterorig$age <=midage_thresh  & masterorig$age > young_thresh),]
# 
# 
# mid_age_sig = matrix (NA, length(49:(dim(master)[2]-1)))
# for (i in 49:(dim(master)[2]-1)) {
#   temp = summary(lm(unlist(master[, i])~master$age + master$age: as.factor(master$geno ) ))
#   temp = temp$coefficients
#   if (temp[3,4]<0.05){
#     mid_age_sig[(i-48)] = "Sig"
#     # print(temp[3,4])
#   }
#   else { mid_age_sig[(i-48)] = "Non_sig" }
# }
# 
# 
# master = masterorig[ (masterorig$age >= midage_thresh ),]
# 
# old_sig = matrix (NA, length(49:(dim(master)[2]-1)))
# for (i in 49:(dim(master)[2]-1)) {
#   temp = summary(lm(unlist(master[, i])~master$age + master$age: as.factor(master$geno ) ))
#   temp = temp$coefficients
#   if (temp[3,4]<0.05){
#     old_sig[(i-48)] = "Sig"
#     # print(temp[3,4])
#   }
#   else { old_sig[(i-48)] = "Non_sig" }
# }
# 
# # now which one is at least significant in one of them
# index_sig = which( old_sig=="Sig" | mid_age_sig =="Sig" |Young_sig=="Sig" ) + 48
# 
# 
# 
# 
# ####### young 
# 
# master = masterorig[masterorig$age <= young_thresh,]
# 
# 
# 
# 
# # for (i in 49:(dim(master)[2]-1)) {
# for (i in index_sig) {
#   temp = summary(lm(unlist(master[, i])~master$age + master$age: as.factor(master$geno ) ))
#   temp = temp$coefficients
#   # if( temp[3,4] <= 0.05){
#     print(i)
#     
#   plt = ggplot(master, aes(x=age, y= unlist(master[, i]) , color=geno, shape=geno)) +
#     geom_point() +
#     geom_smooth(method='lm', fullrange=TRUE, aes(fill = geno)) +
#     theme_linedraw() +
#     theme(text = element_text(size=20))+
#     labs( title = paste0("P-Value = ", temp[3,4]),  x = "Age",   y = paste0( "Mean FA of " , colnames( master)[i]  ) ) 
#   ggsave( paste0("figsmean_young/Mean_",colnames( master)[ i], ".png") )
#   # }  else {print(temp[3,4])}
# }
# 
# 
# ####### midage 
# 
# master = masterorig[ (masterorig$age <=midage_thresh  & masterorig$age > young_thresh),]
# 
# 
# 
# # for (i in 49:(dim(master)[2]-1)) {
# for (i in index_sig) {
#   temp = summary(lm(unlist(master[, i])~master$age + master$age: as.factor(master$geno ) ))
#   temp = temp$coefficients
#   # if( temp[3,4] <= 0.05){
#     print(i)
#     
#     plt = ggplot(master, aes(x=age, y= unlist(master[, i]) , color=geno, shape=geno)) +
#       geom_point() +
#       geom_smooth(method='lm', fullrange=TRUE, aes(fill = geno)) +
#       theme_linedraw() +
#       theme(text = element_text(size=20))+
#       labs( title = paste0("P-Value = ", temp[3,4]),  x = "Age",   y = paste0( "Mean FA of " , colnames( master)[i]  ) ) 
#     ggsave( paste0("figsmean_mid_age/Mean_",colnames( master)[ i], ".png") )
#   # }  else {print(temp[3,4])}
# }
# 
# 
# 
# ####### midage 
# master = masterorig[ (masterorig$age >= midage_thresh ),]
# 
# 
# 
# 
# # for (i in 49:(dim(master)[2]-1)) {
# for (i in index_sig) {
#   temp = summary(lm(unlist(master[, i])~master$age + master$age: as.factor(master$geno ) ))
#   temp = temp$coefficients
#   # if( temp[3,4] <= 0.05){
#     print(i)
#     
#     plt = ggplot(master, aes(x=age, y= unlist(master[, i]) , color=geno, shape=geno)) +
#       geom_point() +
#       geom_smooth(method='lm', fullrange=TRUE, aes(fill = geno)) +
#       theme_linedraw() +
#       theme(text = element_text(size=20))+
#       labs( title = paste0("P-Value = ", temp[3,4]),  x = "Age",   y = paste0( "Mean FA of " , colnames( master)[i]  ) ) 
#     ggsave( paste0("figsmean_old/Mean_",colnames( master)[ i], ".png") )
#   # }  else {print(temp[3,4])}
# }
# 
# 
# 
# 
# 
