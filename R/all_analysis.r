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
					cost=parameters$cost,kernel=as.character(parameters$kernel),gamma=parameters$gamma)
		pred = as.numeric(attr(predict(model, Test[,-dim(Test)[2]],probability = TRUE),"probabilities")[,1])
	}else{
		return("Please select  a ML method!")
	}
	return(pred)
}

fs_order = function(type,Data,Label){
	out = 0
	if(type=="COR"){
		corr = rep(0,(dim(Data)[2]-1))
		for(i in 1:(dim(Data)[2]-1)){
			if((sd(Data[,i])!=0) &&(sd(Label)!=0)){
				corr[i] = cor(Data[,i],Label)
			}
		}
		out = order(abs(corr),decreasing=T)
	}else if(type=="OneR"){
		al_select = c()
		for(numF in 1:(dim(Data)[2]-2)){
			if(numF==1){
				cur_data = Data
				select = OneR(Label ~.,data = optbin(cur_data), verbose = F)
				al_select = c(al_select,which(names(Data)==select$feature))
				cur_data = Data[,-al_select]
			}else{
				cur_data = Data[,-al_select]
			}
			select = OneR(Label ~.,data = optbin(cur_data), verbose = F)
			al_select = c(al_select,which(names(Data)==select$feature))
		}
		out = al_select
	}else{
		cat(paste0("Please type ",paste0("'COR' ",paste0( "or"," 'OneR'!"))),"\n")
	}
	return(out)
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

ncv_perf = matrix(NA,nrow=5,ncol=6)
Rnk = matrix(NA,nrow=dim(Training)[2]-1,ncol=FOLD)
min_numF = c(32, 32, 16, 16, 32)
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
	in_label = rep(-1,dim(Train)[1]);	in_label[which(Train$Label=="R")] = 1
	Rnk[,fold] = fs_order(FS_order,Train,in_label)
	rank = Rnk[,fold]

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
	mx_par = which.max(hyper_grid$auc)
	for(z in 1:(dim(hyper_grid)[2]-1)){
		if(class(hyper_grid[mx_par,z])=="factor"){
			opt_par[fold,z] = as.character(hyper_grid[mx_par,z])
		}else{
			opt_par[fold,z] = hyper_grid[mx_par,z]
		}
	}
	numF = unique(c(seq(min_numF[fold],length(rank),5),length(rank)))
	out_fs_perf = rep(NA,length(numF))
	##### Select features based on outer loop ######
	cat("Select features for testing fold: ",fold, "\n")
	for(nF in 1:length(numF)){
		sel_feature = sort(rank[c(1:numF[nF])])
		out_train = Train[,c(sel_feature,dim(Train)[2])]
		out_test = Test[,c(sel_feature,dim(Test)[2])]
		set.seed(80); rows = sample(nrow(out_train),replace = FALSE)
		out_train <- out_train[rows,];	rm(rows)
		set.seed(50); rows = sample(nrow(out_test),replace = FALSE)
		out_test <- out_test[rows,];	rm(rows)
		pred = ML_method(out_train,out_test,method,hyper_grid[mx_par,])
		test_lab = rep(0,dim(out_test)[1]);	test_lab[which(out_test$Label=="R")] = 1
		out_fs_perf[nF] = as.numeric(auc(test_lab,pred))
		rm(sel_feature);rm(out_train);rm(out_test);rm(pred);rm(test_lab)
	}
	temp_mx = which.max(out_fs_perf)
	test = c(numF[temp_mx]-4,numF[temp_mx]-3,numF[temp_mx]-2,numF[temp_mx]-1,numF[temp_mx],
			numF[temp_mx]+1,numF[temp_mx]+2,numF[temp_mx]+3,numF[temp_mx]+4)
	if(length(which(test>(dim(Training)[2]-1)))!=0){
		test[which(test>(dim(Training)[2]-1))] = dim(Training)[2] - 1
	}
	if(length(which(test<min_numF[fold]))!=0){
		test[which(test<min_numF[fold])] = min_numF[fold]
	}
	test = unique(test)
	Out_fs_perf = matrix(NA,nrow=length(test),ncol=6);	Out_fs_perf[,1] = test
	for(z in 1:length(test)){
		sel_feature = sort(rank[c(1:test[z])])
		out_train = Train[,c(sel_feature,dim(Train)[2])]
		out_test = Test[,c(sel_feature,dim(Test)[2])]
		set.seed(80); rows = sample(nrow(out_train),replace = FALSE)
		out_train <- out_train[rows,];	rm(rows)
		set.seed(50); rows = sample(nrow(out_test),replace = FALSE)
		out_test <- out_test[rows,];	rm(rows)
		pred = ML_method(out_train,out_test,method,hyper_grid[mx_par,])
		test_lab = rep(0,dim(out_test)[1]);	test_lab[which(out_test$Label=="R")] = 1
		Out_fs_perf[z,c(2:6)] = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred,test_lab)
		rm(sel_feature);rm(out_train);rm(out_test);rm(pred);rm(test_lab)
	}

	Mx_id = which.max(Out_fs_perf[,6])[1];	ncv_perf[fold,] = Out_fs_perf[Mx_id,]
	cat("=====================================================","\n")
	cat("FOLD = ",fold,"\n")
	cat("Num Features = ",test[Mx_id],"\n")	
	cat("Performance = ",round(ncv_perf[fold,],digits=4),"\n")
	cat("=====================================================","\n")
	rm(Mx_id);rm(Train);rm(Test);rm(Out_fs_perf);rm(test);rm(hyper_grid)
}
###################################

####### Confirm optimal parameters and features #######
per_mat = matrix(NA,nrow=FOLD,ncol=10)
for(par in 1:FOLD){
	temp_perf = matrix(NA,nrow=5,ncol=5)
	sel_feature = sort(Rnk[c(1:ncv_perf[par,1]),par])
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
		Train <- Train[,c(sel_feature,dim(Train)[2])]
		Test <- Test[,c(sel_feature,dim(Test)[2])]
		pred = ML_method(Train,Test,method,opt_par[par,])
		test_lab = rep(0,dim(Test)[1]);	test_lab[which(Test$Label=="R")] = 1
		temp_perf[fold,] = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred,test_lab)
		rm(pred);rm(Train);rm(Test);rm(test_lab)
		cat("=====================================================","\n")
		cat("FOLD = ",fold,"\n")
		cat("Performance = ",round(ncv_perf[fold,],digits=4),"\n")
		cat("=====================================================","\n")
	}#end_5-fold
	per_mat[par,] = c(mean(temp_perf[,1]),mean(temp_perf[,2]),mean(temp_perf[,3]),mean(temp_perf[,4]),mean(temp_perf[,5]),
					sd(temp_perf[,1]),sd(temp_perf[,2]),sd(temp_perf[,3]),sd(temp_perf[,4]),sd(temp_perf[,5]))
	rm(temp_perf);rm(sel_feature)
}#end_par
###################################
opt_idx = which.max(per_mat[,5])
opt_numF = ncv_perf[opt_idx,1]
Opt_par = opt_par[opt_idx,]
fs_perf = per_mat[opt_idx,]
sel_feature = sort(Rnk[c(1:ncv_perf[opt_idx,1]),opt_idx])
final_perf = per_mat[opt_idx,]

####### Independent testing #######
indep_Train <- Training
indep_Test <- Testing
set.seed(80); rows = sample(nrow(indep_Train),replace = FALSE)
indep_Train <- indep_Train[rows,];	rm(rows)
set.seed(50); rows = sample(nrow(indep_Test),replace = FALSE)
indep_Test <- indep_Test[rows,];	rm(rows)
test_lab = rep(0,dim(indep_Test)[1]);	test_lab[which(indep_Test$Label=="R")] = 1
pred0 = ML_method(indep_Train,indep_Test,method,Opt_par)
indep_perf_NO_FS = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred0,test_lab)
rm(indep_Train);rm(indep_Test)
set.seed(80); rows = sample(nrow(Training),replace = FALSE)
indep_Train <- Training[rows,c(sel_feature,dim(Training)[2])]
set.seed(50); rows = sample(nrow(Testing),replace = FALSE)
indep_Test <- Testing[rows,c(sel_feature,dim(Testing)[2])]
pred1 = ML_method(indep_Train,indep_Test,method,Opt_par)
indep_perf_FS = eval_fun(length(which(test_lab==1)),length(which(test_lab==0)),pred1,test_lab)

###################################
cat(" ","\n",
	"*************************************","\n",
	"Performance of 5-fold cross validation with feature selection:","\n",
	"Sensitivity: ",round(fs_perf[1],digits=4),"±",round(fs_perf[6],digits=4),"\n",
	"Specificity: ",round(fs_perf[2],digits=4),"±",round(fs_perf[7],digits=4),"\n",
	"Accuracy: ",round(fs_perf[3],digits=4),"±",round(fs_perf[8],digits=4),"\n",
	"Matthews correlation coefficient: ",round(fs_perf[4],digits=4),"±",round(fs_perf[9],digits=4),"\n",
	"AUC: ",round(fs_perf[5],digits=4),"±",round(fs_perf[10],digits=4),"\n",
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
