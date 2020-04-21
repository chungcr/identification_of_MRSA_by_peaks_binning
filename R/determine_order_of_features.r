rm(list=ls(all=TRUE))
library(OneR)

output_path = paste0(getwd(),"/")

data = read.csv("https://raw.githubusercontent.com/chungcr/identification_of_MRSA_by_peaks_binning/master/data/Training%20data.csv")
label = rep(-1,dim(data)[1])
label[which(data$Label=="R")] = 1


########## Order by Pearson's correlation coefficients ##########
corr = rep(NA,(dim(data)[2]-1))
for(i in 1:(dim(data)[2]-1)){
	corr[i] = cor(data[,i],label)
}
COR = order(abs(corr),decreasing=T)

########## Order by OneR ##########
al_select = c()
for(numF in 1:(dim(data)[2]-2)){
	if(numF==1){
		cur_data = data
		select = OneR(Label ~.,data = optbin(cur_data), verbose = F)
		al_select = c(al_select,which(names(data)==select$feature))
		cur_data = data[,-al_select]
	}else{
		cur_data = data[,-al_select]
	}
	select = OneR(Label ~.,data = optbin(cur_data), verbose = F)
	al_select = c(al_select,which(names(data)==select$feature))
}
OneR = al_select

final = data.frame(COR,OneR)
write.csv(final,paste0(output_path,"FS_order.csv"),row.names=F)