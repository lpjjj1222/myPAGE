setwd("/Users/lpjjj/Downloads")

filtered<-read.table('/Users/lpjjj/Downloads/filtered_table.txt')
#filtered<-read.table('/Users/lpjjj/Downloads/scPAGE-main-supplement/model2/for_training/sepsis_189400.series.matrix.txt')
#filtered<-t(filtered)
label<-read.table('/Users/lpjjj/Downloads/scPAGE-main-supplement/model2/label/label_189400.txt')
label<-label[-1,]
genes<-colnames(filtered)
genes

result_df <- data.frame(
  'A' = character(0),
  'b' = character(0),
  'A大于Bin_N' = numeric(0),
  'A小于Bin_S' = numeric(0),
  stringsAsFactors = FALSE)


for (i in 1:(ncol(filtered) - 1)) {
  for (j in (i + 1):ncol(filtered)) {
    gene_a<-genes[i]
    gene_b<-genes[j]
    print('************')
    print(gene_a)
    print(gene_b)
    #新的基因对初始化
    positive_count<-0
    negative_count<-0
    for (sample in 1:nrow(filtered)) {
      labels =label[sample]
      A <- filtered[sample,i]
      B <- filtered[sample,j]

      # 计算基因差值
      diff_values <- A - B
      # 统计差值在不同标签下的样本数量
      if(diff_values > 0 & labels == 1) {
        positive_count<-positive_count+1
      }else if (diff_values < 0 & labels == 0) {
        negative_count<-negative_count+1
      }
      
      # 将统计结果添加到结果数据框

    }
    print(positive_count)
    print(negative_count)
    result_df <- rbind(result_df, data.frame('A' = gene_a,'B' = gene_b, 'A大于Bin_N' = positive_count,'A小于Bin_S' = negative_count))
}
  
}

#-------------------------------------------
saveRDS(result_df,file = 'result_df.rds')





#有多少个S样本（SEPSIS）
S<-length(which(label == 1))
#有多少个N样本(NORMAL)
N<-length(which(label==0))

data <- result_df
colnames(result_df)<-c('A','B','A大于BinS','A小于BinN')

#load("data.Rdata")
#fisher.t.test#####



fishfun=function(x){
  x1=as.numeric(x[1])
  x2=as.numeric(x[2])
  d2=matrix(c(N-x2,x1,x2,S-x1),ncol=2,nrow=2)
  t.res=fisher.test(d2)$p.value
}
ratio_label=apply(data[,c(3,4)],1,fishfun )


index_9=which(ratio_label<10^(-9)) #434个
rev_pair_9=data[index,]
index_8=which(ratio_label<10^(-8)) #784个
rev_pair_8=data[index,]
index_7=which(ratio_label<10^(-7))#1387个
rev_pair_7=data[index,]

rev_pair_9<-rev_pair_9[,c(1,2)]
rev_pair_8<-rev_pair_8[,c(1,2)]
rev_pair_7<-rev_pair_7[,c(1,2)]








