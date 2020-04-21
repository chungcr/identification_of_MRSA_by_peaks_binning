########## Import packages ##########
library(caret)
library(pROC)
library(OptimalCutpoints)
library(rpart)	#DT
library(ranger)	#RF
library(kknn)	#KNN
library(e1071)	#SVM

rm(list=ls(all=TRUE))

########## Set data path ##########
FOLD = 5
method = "DT"		#Choose ML method ("DT","RF","KNN","SVM")
FS_order = "COR"		#Choose order of feature for forward feature selection ("COR","OneR")

Training = read.csv("https://raw.githubusercontent.com/chungcr/identification_of_MRSA_by_peaks_binning/master/data/Training%20data.csv")
Testing = read.csv("https://raw.githubusercontent.com/chungcr/identification_of_MRSA_by_peaks_binning/master/data/Testing%20data.csv")

ML_grid = function(ML){
	if(ML == "DT"){
		hyper_grid <- expand.grid(
			minsplit = seq(5,40,10),
			minbucket = round(seq(5,40,5)/3),
			cp = c(0.01, 0.1, 0.5),
			maxdepth = seq(5,30,5),
			auc = 0
		)
	}else if(ML == "RF"){
		hyper_grid <- expand.grid(
			mtry = c(16,32,128),
			ntre = c(100,500,1000,1500),
			node_size = c(1,5,15),
			sample_size = c(0.5, 0.632,1),
			auc = 0
		)
	}else if(ML == "KNN"){
		hyper_grid <- expand.grid(
			k = seq(2,70,5),
			kernel = c("rectangular","triangular","gaussian","optimal"),
			auc = 0
		)
	}else if(ML == "SVM"){
		hyper_grid <- expand.grid(
			cost = c(0.1, 1, 10, 100, 1000),
			gamma = c(1/534,0.5, 1, 2),
			kernel = c("linear","polynomial","radial","sigmoid"),
			auc = 0
		)
	}else{
		return("Please select  a ML method!")
	}
	return(hyper_grid)
}


ML_method = function(Train,Test,ML,parameters){
	if(ML == "DT"){
		set.seed(10);model = rpart(as.factor(Label)~.,data = Train,method="class",
				minsplit=parameters$minsplit,minbucket=parameters$minbucket,
				cp=parameters$cp, maxdepth=parameters$maxdepth)
		pred = as.numeric(predict(model,newdata = Test[,-dim(Test)[2]],type="prob")[,1])
	}else if(ML == "RF"){
		set.seed(10);model = ranger(Label~.,data = Train,
					num.trees = parameters$ntre,
					mtry = parameters$mtry,
					min.node.size = parameters$node_size,
					sample.fraction = parameters$sample_size,probability = TRUE,seed=10)
		pred = predict(model,data = Test[,-dim(Test)[2]])$predictions[,1]
	}else if(ML == "KNN"){
		set.seed(10);model = kknn(Train$Label~.,Train[,-dim(Train)[2]],Test[,-dim(Test)[2]], 
					k = parameters$k,kernel = as.character(parameters$kernel))
		pred = model$prob[,1]
	}else if(ML == "SVM"){
		set.seed(10);model = svm(Train[,-dim(Train)[2]],as.factor(Train$Label),probability = TRUE,
					cost=parameters$cost,kernel=as.character(parameters$kernel),gamma=para$gamma)
		pred = as.numeric(attr(predict(model, Test[,-dim(Test)[2]],probability = TRUE),"probabilities")[,1])
	}else{
		return("Please select  a ML method!")
	}
	return(pred)
}

eval_fun = function(n_pos,n_neg,pred,true){
	Data = data.frame(pred,true)
	optcut = optimal.cutpoints(X = "pred",status = "true" ,methods = "SpEqualSe",data=Data, tag.healthy = 0)
	fp = as.numeric(optcut$SpEqualSe$Global$optimal.cutoff$FP)[length(optcut$SpEqualSe$Global$optimal.cutoff$`cutoff`)]
	fn = as.numeric(optcut$SpEqualSe$Global$optimal.cutoff$FN)[length(optcut$SpEqualSe$Global$optimal.cutoff$`cutoff`)]
	tp = as.numeric(n_pos - fn);			tn = as.numeric(n_neg - fp)
	cv_sen = tp/n_pos;				cv_spe = tn/n_neg
	cv_acc = (tp+tn)/(n_pos+n_neg);		cv_mcc = (tp*tn-fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
	cv_auc = as.numeric(auc(true,pred))
	c(cv_sen,cv_spe,cv_acc,cv_mcc,cv_auc)
}	

tic = Sys.time()

label = rep(NA,dim(Training)[1]);	label[which(Training$Label=="R")] = 1;	label[which(Training$Label=="S")] = 0
pos = c(1:length(which(Training$Label=="R")));	neg = c(1:length(which(Training$Label=="S")))
n = 10;	set.seed(n);	ffold_pos = createFolds(pos,5);	set.seed(n+200);	ffold_neg = createFolds(neg,5)

para_name = names(ML_grid(method));	para_name = para_name[-length(para_name)]
opt_par = matrix(NA,nrow=5,ncol=length(para_name))
opt_par = data.frame(opt_par); names(opt_par) = para_name
ncv_perf = matrix(NA,nrow=5,ncol=5)

for(fold in 1:5){
	cat("====================","\n")
	cat("START FOLD = ",fold,"\n")
	train_pos = c();	train_neg = c()
	test_pos = ffold_pos[[fold]];	test_neg = ffold_neg[[fold]]
	temp_fold = c(1:FOLD);	temp_fold = temp_fold[-fold]
	for(jj in 1:length(temp_fold)){
		train_pos = c(train_pos,ffold_pos[[temp_fold[jj]]]);
		train_neg = c(train_neg,ffold_neg[[temp_fold[jj]]]);
	}
	Train = Training[c(train_pos,length(which(Training$Label=="R"))+train_neg),]
	Test = Training[c(test_pos,length(which(Training$Label=="R"))+test_neg),]
	sub_label = rep(NA,dim(Train)[1]);	sub_label[which(Train$Label=="R")] = 1;	sub_label[which(Train$Label=="S")] = 0
	sub_pos = c(1:length(which(Train$Label=="R")));	sub_neg = c(1:length(which(Train$Label=="S")))
	n = 10;	set.seed(n);	sub_ffold_pos = createFolds(sub_pos,FOLD);
			set.seed(n+200);	sub_ffold_neg = createFolds(sub_neg,FOLD)
	hyper_grid = ML_grid(method)
	for(par in 1:nrow(hyper_grid)){
		temp_perf = rep(NA,FOLD)
		for(sub_fold in 1:FOLD){
			cat("subFOLD = ",sub_fold ,"\n")
			sub_train_pos = c();	sub_train_neg = c()
			sub_test_pos = sub_ffold_pos[[sub_fold]];	sub_test_neg = sub_ffold_neg[[sub_fold]]
			temp_fold = c(1:FOLD);	temp_fold = temp_fold[-sub_fold]
			for(jj in 1:length(temp_fold)){
				sub_train_pos = c(sub_train_pos,sub_ffold_pos[[temp_fold[jj]]]);
				sub_train_neg = c(sub_train_neg,sub_ffold_neg[[temp_fold[jj]]]);
			}
			sub_Train = Train[c(sub_train_pos,length(which(Train$Label=="R"))+sub_train_neg),]
			sub_Test = Train[c(sub_test_pos,length(which(Train$Label=="R"))+sub_test_neg),]
			set.seed(80); rows = sample(nrow(sub_Train),replace = FALSE)
			sub_Train <- sub_Train[rows,];	rm(rows)
			set.seed(50); rows = sample(nrow(sub_Test),replace = FALSE)
			sub_Test <- sub_Test[rows,];	rm(rows)
			pred = ML_method(sub_Train,sub_Test,method,hyper_grid[par,])
			test_lab = rep(0,dim(sub_Test)[1]);	test_lab[which(sub_Test$Label=="R")] = 1
				temp_perf[sub_fold] = as.numeric(auc(test_lab,pred))
				cat("subfold_AUC = ",temp_perf[sub_fold] ,"\n")
				rm(pred);rm(test_lab);rm(sub_Train);rm(sub_Test)
		}
		hyper_grid$auc[par] = mean(temp_perf)
		cat("*******************","\n")
		cat("FOLD = ",fold,"\n")
		cat("Finished = ",par,"/",dim(hyper_grid)[1],"\n")
		cat("AUC = ",hyper_grid$auc[par],"\n");	rm(temp_perf)
		cat("*******************","\n")
	}#end_parameter
	set.seed(80); rows = sample(nrow(Train),replace = FALSE)
	Train <- Train[rows,];	rm(rows)
	set.seed(50); rows = sample(nrow(Test),replace = FALSE)
	Test <- Test[rows,];	rm(rows)
	mx_id = which.max(hyper_grid$auc)
	for(z in 1:(dim(hyper_grid)[2]-1)){
		if(class(hyper_grid[mx_id,z])=="factor"){
			opt_par[fold,z] = as.character(hyper_grid[mx_id,z])
		}else{
			opt_par[fold,z] = hyper_grid[mx_id,z]
		}
	}
	pred = ML_method(Train,Test,method,hyper_grid[mx_id,])
	test_lab = rep(0,dim(Test)[1]);	test_lab[which(Test$Label=="R")] = 1
	ncv_perf[fold,] = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred,test_lab)
	cat("End FOLD = ",fold,"\n")
	cat("====================","\n")
	rm(pred);rm(Train);rm(Test);rm(test_lab);rm(hyper_grid)
}
###################################
uni_par = unique.matrix(opt_par);	rm(opt_par)

####### Confirm optimal parameters #######
per_mat = matrix(NA,nrow=dim(uni_par)[1],ncol=10)
for(par in 1:dim(uni_par)[1]){
	ncv_perf = matrix(NA,nrow=5,ncol=5)
	for(fold in 1:FOLD){
		train_pos = c();	train_neg = c()
		test_pos = ffold_pos[[fold]];	test_neg = ffold_neg[[fold]]
		temp_fold = c(1:FOLD);	temp_fold = temp_fold[-fold]
		for(jj in 1:length(temp_fold)){
			train_pos = c(train_pos,ffold_pos[[temp_fold[jj]]]);
			train_neg = c(train_neg,ffold_neg[[temp_fold[jj]]]);
		}
		Train = Training[c(train_pos,length(which(Training$Label=="R"))+train_neg),]
		Test = Training[c(test_pos,length(which(Training$Label=="R"))+test_neg),]
		set.seed(80); rows = sample(nrow(Train),replace = FALSE)
		Train <- Train[rows,];	rm(rows)
		set.seed(50); rows = sample(nrow(Test),replace = FALSE)
		Test <- Test[rows,];	rm(rows)
		pred = ML_method(Train,Test,method,uni_par[par,])
		test_lab = rep(0,dim(Test)[1]);	test_lab[which(Test$Label=="R")] = 1
		ncv_perf[fold,] = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred,test_lab)
		rm(pred);rm(Train);rm(Test);rm(test_lab)
		cat("=====================================================","\n")
		cat("FOLD = ",fold,"\n")
		cat("Performance = ",round(ncv_perf[fold,],digits=4),"\n")
		cat("=====================================================","\n")
	
	}#end_5-fold
	per_mat[par,] = c(mean(ncv_perf[,1]),mean(ncv_perf[,2]),mean(ncv_perf[,3]),mean(ncv_perf[,4]),mean(ncv_perf[,5]),
					sd(ncv_perf[,1]),sd(ncv_perf[,2]),sd(ncv_perf[,3]),sd(ncv_perf[,4]),sd(ncv_perf[,5]))
	rm(ncv_perf)
}#end_par
###################################
mx_id = which.max(per_mat[,5])
final_perf = per_mat[mx_id,]
opt_par = uni_par[mx_id,]
opt_ML_par = opt_par

####### Start feature selection #######
fs_order = read.csv("https://raw.githubusercontent.com/chungcr/identification_of_MRSA_by_peaks_binning/master/data/FS_order.csv")
if(FS_order == "OneR"){
	rank = fs_order$OneR
}else{
	rank = fs_order$COR
}
ncv_perf = matrix(NA,nrow=5,ncol=6)
for(fold in 1:5){
	cat("====================","\n")
	cat("START FOLD = ",fold,"\n")
	train_pos = c();	train_neg = c()
	test_pos = ffold_pos[[fold]];	test_neg = ffold_neg[[fold]]
	temp_fold = c(1:FOLD);	temp_fold = temp_fold[-fold]
	for(jj in 1:length(temp_fold)){
		train_pos = c(train_pos,ffold_pos[[temp_fold[jj]]]);
		train_neg = c(train_neg,ffold_neg[[temp_fold[jj]]]);
	}
	Train = Training[c(train_pos,length(which(Training$Label=="R"))+train_neg),]
	Test = Training[c(test_pos,length(which(Training$Label=="R"))+test_neg),]
	sub_label = rep(NA,dim(Train)[1]);	sub_label[which(Train$Label=="R")] = 1;	sub_label[which(Train$Label=="S")] = 0
	sub_pos = c(1:length(which(Train$Label=="R")));	sub_neg = c(1:length(which(Train$Label=="S")))
	n = 10;	set.seed(n);	sub_ffold_pos = createFolds(sub_pos,FOLD);
			set.seed(n+200);	sub_ffold_neg = createFolds(sub_neg,FOLD)

	if(method == "RF"){
		min_mtry = opt_par$mtry
		try = ceiling((length(rank)-min_mtry)/5)
	}else{
		try = round(length(rank)/5)
	}
	temp_perf = matrix(NA,nrow=try,ncol=FOLD)
	for(sub_fold in 1:FOLD){
		sub_train_pos = c();	sub_train_neg = c()
		sub_test_pos = sub_ffold_pos[[sub_fold]];	sub_test_neg = sub_ffold_neg[[sub_fold]]
		temp_fold = c(1:FOLD);	temp_fold = temp_fold[-sub_fold]
		for(jj in 1:length(temp_fold)){
			sub_train_pos = c(sub_train_pos,sub_ffold_pos[[temp_fold[jj]]]);
			sub_train_neg = c(sub_train_neg,sub_ffold_neg[[temp_fold[jj]]]);
		}
		Sub_Train = Train[c(sub_train_pos,length(which(Train$Label=="R"))+sub_train_neg),]
		Sub_Test = Train[c(sub_test_pos,length(which(Train$Label=="R"))+sub_test_neg),]
		set.seed(80); rows = sample(nrow(Sub_Train),replace = FALSE)
		Sub_Train <- Sub_Train[rows,];	rm(rows)
		set.seed(50); rows = sample(nrow(Sub_Test),replace = FALSE)
		Sub_Test <- Sub_Test[rows,];	rm(rows)
		for(par in 1:try){
			if(method == "RF"){
				if(par==1){
					sel_feature = rank[c(1:(min_mtry))][which(!is.na(rank[c(1:(min_mtry))]))]
					sub_Train = Sub_Train[,c(sel_feature,dim(Sub_Train)[2])]
					sub_Test = Sub_Test[,c(sel_feature,dim(Sub_Test)[2])]
				}else{
					sel_feature = rank[c(1:(min_mtry+5*(par-1)))][which(!is.na(rank[c(1:(min_mtry+5*(par-1)))]))]
					sub_Train = Sub_Train[,c(sel_feature,dim(Sub_Train)[2])]
					sub_Test = Sub_Test[,c(sel_feature,dim(Sub_Test)[2])]
				}
				pred = ML_method(sub_Train,sub_Test,method,opt_par)
				test_lab = rep(0,dim(sub_Test)[1]);	test_lab[which(sub_Test$Label=="R")] = 1
				temp_perf[par,sub_fold] = as.numeric(auc(test_lab,pred))
				cat("*********************************************","\n")
				cat("FOLD = ",fold,"\n")
				cat("subFOLD = ",sub_fold ,"\n")
				cat(par,"/",try,"\n")
				cat("subfold_AUC = ",temp_perf[par,sub_fold] ,"\n")
				cat("*********************************************","\n")
				rm(pred);rm(test_lab);rm(sub_Train);rm(sub_Test);rm(sel_feature)
			}else{
				sel_feature = rank[c(1:(par*5))][which(!is.na(rank[c(1:(par*5))]))]
				sub_Train = Sub_Train[,c(sel_feature,dim(Sub_Train)[2])]
				sub_Test = Sub_Test[,c(sel_feature,dim(Sub_Test)[2])]
				pred = ML_method(sub_Train,sub_Test,method,opt_par)
				test_lab = rep(0,dim(sub_Test)[1]);	test_lab[which(sub_Test$Label=="R")] = 1
				temp_perf[par,sub_fold] = as.numeric(auc(test_lab,pred))
				cat("*********************************************","\n")
				cat("FOLD = ",fold,"\n")
				cat("subFOLD = ",sub_fold ,"\n")
				cat(par,"/",try,"\n")
				cat("subfold_AUC = ",temp_perf[par,sub_fold] ,"\n")
				cat("*********************************************","\n")
				rm(pred);rm(test_lab);rm(sub_Train);rm(sub_Test);rm(sel_feature)
			}	
		}
	}#end_parameter
	eva = rep(0,try)
	for(j in 1:try){
		eva[j] = mean(temp_perf[j,], na.rm = TRUE)
	}
	set.seed(80); rows = sample(nrow(Train),replace = FALSE)
	Train <- Train[rows,];	rm(rows)
	set.seed(50); rows = sample(nrow(Test),replace = FALSE)
	Test <- Test[rows,];	rm(rows)
	mx_id = which.max(eva)[1]	;rm(eva)
	test = c(5*mx_id-3,5*mx_id-2,5*mx_id-1,5*mx_id,5*mx_id+1,5*mx_id+2)
	Temp_perf = matrix(NA,nrow=length(test),ncol=5)
	for(z in 1:length(test)){
		Ttrain <- Train[,c(rank[c(1:test[z])],dim(Sub_Train)[2])]
		Ttest <- Test[,c(rank[c(1:test[z])],dim(Sub_Train)[2])]
		pred = ML_method(Ttrain,Ttest,method,opt_par)
		test_lab = rep(0,dim(Ttest)[1]);	test_lab[which(Ttest$Label=="R")] = 1
		Temp_perf[z,] = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred,test_lab)
		rm(pred);rm(Ttrain);rm(Ttest);rm(test_lab)
	}
	Mx_id = which.max(Temp_perf[,5])[1]
	ncv_perf[fold,] = c(test[Mx_id],Temp_perf[Mx_id,])
	cat("=====================================================","\n")
	cat("FOLD = ",fold,"\n")
	cat("Num Features = ",test[Mx_id],"\n")	
	cat("Performance = ",round(ncv_perf[fold,],digits=4),"\n")
	cat("=====================================================","\n")
	rm(mx_id);rm(Temp_perf);rm(Mx_id);rm(test);rm(Train);rm(Test)
}
###################################
Uni_par = unique(ncv_perf[,1]);

####### Confirm optimal features #######
per_mat = matrix(NA,nrow=length(Uni_par),ncol=10)
for(par in 1:length(Uni_par)){
	ncv_perf = matrix(NA,nrow=5,ncol=5)
	for(fold in 1:FOLD){
		train_pos = c();	train_neg = c()
		test_pos = ffold_pos[[fold]];	test_neg = ffold_neg[[fold]]
		temp_fold = c(1:FOLD);	temp_fold = temp_fold[-fold]
		for(jj in 1:length(temp_fold)){
			train_pos = c(train_pos,ffold_pos[[temp_fold[jj]]]);
			train_neg = c(train_neg,ffold_neg[[temp_fold[jj]]]);
		}
		Train = Training[c(train_pos,length(which(Training$Label=="R"))+train_neg),]
		Test = Training[c(test_pos,length(which(Training$Label=="R"))+test_neg),]
		set.seed(80); rows = sample(nrow(Train),replace = FALSE)
		Train <- Train[rows,];	rm(rows)
		set.seed(50); rows = sample(nrow(Test),replace = FALSE)
		Test <- Test[rows,];	rm(rows)
		Train <- Train[,c(rank[c(1:Uni_par[par])],dim(Train)[2])]
		Test <- Test[,c(rank[c(1:Uni_par[par])],dim(Test)[2])]
		pred = ML_method(Train,Test,method,opt_par)
		test_lab = rep(0,dim(Test)[1]);	test_lab[which(Test$Label=="R")] = 1
		ncv_perf[fold,] = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred,test_lab)
		rm(pred);rm(Train);rm(Test);rm(test_lab)
		cat("=====================================================","\n")
		cat("FOLD = ",fold,"\n")
		cat("Performance = ",round(ncv_perf[fold,],digits=4),"\n")
		cat("=====================================================","\n")
	}#end_5-fold
	per_mat[par,] = c(mean(ncv_perf[,1]),mean(ncv_perf[,2]),mean(ncv_perf[,3]),mean(ncv_perf[,4]),mean(ncv_perf[,5]),
					sd(ncv_perf[,1]),sd(ncv_perf[,2]),sd(ncv_perf[,3]),sd(ncv_perf[,4]),sd(ncv_perf[,5]))
	rm(ncv_perf)
}#end_par
###################################
opt_numF = Uni_par[which.max(per_mat[,5])]
fs_perf = per_mat[which.max(per_mat[,5]),]

####### Independent testing #######
indep_Train <- Training
indep_Test <- Testing
set.seed(80); rows = sample(nrow(indep_Train),replace = FALSE)
indep_Train <- indep_Train[rows,];	rm(rows)
set.seed(50); rows = sample(nrow(indep_Test),replace = FALSE)
indep_Test <- indep_Test[rows,];	rm(rows)
test_lab = rep(0,dim(indep_Test)[1]);	test_lab[which(indep_Test$Label=="R")] = 1
pred0 = ML_method(indep_Train,indep_Test,method,opt_par)
indep_perf_NO_FS = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred0,test_lab)
rm(indep_Train);rm(indep_Test)
set.seed(80); rows = sample(nrow(Training),replace = FALSE)
indep_Train <- Training[rows,c(rank[c(1:opt_numF)],dim(Training)[2])]
set.seed(50); rows = sample(nrow(Testing),replace = FALSE)
indep_Test <- Testing[rows,c(rank[c(1:opt_numF)],dim(Testing)[2])]
pred1 = ML_method(indep_Train,indep_Test,method,opt_par)
indep_perf_FS = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred1,test_lab)

###################################

cat(" ","\n",
	"*************************************","\n",
	"Performance of 5-fold cross validation without feature selection:","\n",
	"Sensitivity: ",round(final_perf[1],digits=4),"กำ",round(final_perf[6],digits=4),"\n",
	"Specificity: ",round(final_perf[2],digits=4),"กำ",round(final_perf[7],digits=4),"\n",
	"Accuracy: ",round(final_perf[3],digits=4),"กำ",round(final_perf[8],digits=4),"\n",
	"Matthews correlation coefficient: ",round(final_perf[4],digits=4),"กำ",round(final_perf[9],digits=4),"\n",
	"AUC: ",round(final_perf[5],digits=4),"กำ",round(final_perf[10],digits=4),"\n",
	"*************************************","\n",
	"Performance of 5-fold cross validation with feature selection:","\n",
	"Sensitivity: ",round(fs_perf[1],digits=4),"กำ",round(fs_perf[6],digits=4),"\n",
	"Specificity: ",round(fs_perf[2],digits=4),"กำ",round(fs_perf[7],digits=4),"\n",
	"Accuracy: ",round(fs_perf[3],digits=4),"กำ",round(fs_perf[8],digits=4),"\n",
	"Matthews correlation coefficient: ",round(fs_perf[4],digits=4),"กำ",round(fs_perf[9],digits=4),"\n",
	"AUC: ",round(fs_perf[5],digits=4),"กำ",round(fs_perf[10],digits=4),"\n",
	"*************************************","\n",
	"Performance of independent testing without feature selection:","\n",
	"Number of features: ",dim(Training)[2]-1,"\n",
	"Sensitivity: ",round(indep_perf_NO_FS[1],digits=4),"\n",
	"Specificity: ",round(indep_perf_NO_FS[2],digits=4),"\n",
	"Accuracy: ",round(indep_perf_NO_FS[3],digits=4),"\n",
	"Matthews correlation coefficient: ",round(indep_perf_NO_FS[4],digits=4),"\n",
	"AUC: ",round(indep_perf_NO_FS[5],digits=4),"\n",
	"*************************************","\n",
	"Performance of independent testing with feature selection:","\n",
	"Number of features: ",dim(indep_Train)[2]-1,"\n",
	"Sensitivity: ",round(indep_perf_FS[1],digits=4),"\n",
	"Specificity: ",round(indep_perf_FS[2],digits=4),"\n",
	"Accuracy: ",round(indep_perf_FS[3],digits=4),"\n",
	"Matthews correlation coefficient: ",round(indep_perf_FS[4],digits=4),"\n",
	"AUC: ",round(indep_perf_FS[5],digits=4),"\n",
	"*************************************","\n")

toc = Sys.time()
toc - tic

