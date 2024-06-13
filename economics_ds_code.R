######################
#Project code for the 
#Economics for Data Science Exam
#
#Pesaresi Nicola
#m. 842468
######################

dt<- read.csv("marketing_campaign.csv", sep=";", header=T)
dt<-dt[-2241,]
#Explorative analysis and preprocessing
library(ggplot2)

str(dt)
summary(dt)

dt$Response<- as.integer(dt$Response)

dt$Dt_Customer<- round(difftime(Sys.Date(), dt$Dt_Customer, units = "days"),0)
dt$Dt_Customer<- as.integer(dt$Dt_Customer)
ggplot() +
  geom_histogram(aes(Year_Birth), data=dt)
dt<- dt[- which(dt$Year_Birth<1920),] #remove wrong values
dt<- dt[- which(dt$Marital_Status %in% c("Absurd","Alone","YOLO")),]

dt<- dt[- which.max(dt$Income),] #remove wrong values

sum(is.na(dt))
colSums(is.na(dt))
#dt<- na.omit(dt) #remove NAs

#input income NAs
na.rows<- which(rowSums(is.na(dt)) > 0)
income.fit<- lm(Income ~ ., dt[-na.rows,])
summary(income.fit)
income.fit.step<- lm(Income ~ Education + Teenhome + Kidhome,
                     dt[-na.rows,])
summary(income.fit.step)
incms<-predict(income.fit.step, 
               dt[na.rows, c("Education","Teenhome","Kidhome")])
incms<- round(incms,0)
dt[na.rows,"Income"]<- incms
str(dt)



ggplot() +
  geom_histogram(aes(Income), data=dt)


dt.num<- dt[, sapply(dt, is.numeric)]
par(mar = c(11,3,2,2))
boxplot(scale(dt.num)[,c(2:17)], main="", las=2)

#GGally::ggpairs(dt.num)
ggcorrplot::ggcorrplot(cor(dt.num[,-c(1,25,26)]), lab=F, type="upper")

table(dt$Response)

cmp_sold<- data.frame(Campaign=colnames(dt)[c(24,25,21,22,23,29)],
                      Sold=colSums(dt[,c(24,25,21,22,23,29)]),
                      ncol=2)

ggplot(data=cmp_sold) +
  geom_bar(aes(Campaign, weight=Sold, fill="red")) +
  guides(fill="none") +
  theme_bw()

prod_sold<- data.frame(Prod=colnames(dt)[10:15],
                      Sold=colSums(dt[,10:15]),
                      ncol=2)

ggplot(data=prod_sold) +
  geom_bar(aes(Prod, weight=Sold))

where_sold<- data.frame(Where=colnames(dt)[16:19],
                       Sold=colSums(dt[,16:19]),
                       ncol=2)

ggplot(data=where_sold) +
  geom_bar(aes(Where, weight=Sold), fill="cornflowerblue") +
  guides(fill="none") + 
  theme_bw()

#######
#Study of the causal impact of 
#the variable Kidhome on the response
library(xtable)

ggplot(data=dt) +
  geom_point(aes(x=Kidhome, y=Response))

Kidhome.bi<- dt$Kidhome
Kidhome.bi[which(dt$Kidhome >0)]<- 1

table(dt$Kidhome, dt$Response)
table(Kidhome.bi, dt$Response)

cor(dt$Kidhome, dt$Response)
cor(Kidhome.bi, dt$Response)

mean(dt[which(dt$Kidhome==0), "Response"]) #0,17
mean(dt[which(dt$Kidhome==1), "Response"]) #0.12
mean(dt[which(dt$Kidhome==2), "Response"]) #0.04

mean(dt[which(Kidhome.bi==0), "Response"]) #0,17
mean(dt[which(Kidhome.bi==1), "Response"]) #0.12

t.test(dt$Response ~ Kidhome.bi > 0, alternative="greater")
#refuted equal means hyphothesis

mod1<- lm(dt$Response ~ dt$Kidhome)
mod2<- lm(dt$Response ~ Kidhome.bi)
summary(mod1)
summary(mod2)

mod3<- lm(Response ~ Year_Birth + Dt_Customer +
            Education + Kidhome , data=dt)
summary(mod3)

mod4<- lm(Response ~ ., data=dt)
summary(mod4)

library(MatchIt)
dt.match<- dt
dt.match$Kidhome<- Kidhome.bi
match<- matchit(Kidhome ~ ., data=dt.match, method="nearest")
matched.data<- match.data(match)

cobalt::love.plot(match, binary="std")

mod.matching<- lm(Response ~ Kidhome, data=matched.data)
summary(mod.matching)

dt.match$propensity_score <- match$distance
ggplot(dt.match, aes(x = propensity_score, fill = factor(Kidhome))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution Before Matching",
       fill = "Treatment") +
  theme_minimal()

# Propensity scores after matching
matched.data$propensity_score <- match$distance[match$weights > 0]
ggplot(matched.data, aes(x = propensity_score, fill = factor(Kidhome))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution After Matching",
       fill = "Treatment") +
  theme_minimal()


#######
#Market segmentation using Self Organizion Maps
#with the kohonen package

library(kohonen)
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6,
                                                  alpha=alpha)[n:1]}

dt.dummies<- fastDummies::dummy_cols(dt, remove_first_dummy = T,
                                     remove_selected_columns = T)
dt.subs<- dt[,c(10:15,21:25,29)] #subset to be used for clustering
dt.scaled<- as.matrix(scale(dt.subs))

5*sqrt(nrow(dt)) #number of neurons. ~15*15
som.grid <- somgrid(xdim = 15, ydim = 15, topo = "hexagonal")

som.model <- som(dt.scaled,
                 grid = som.grid, 
                 rlen = 10000,
                 alpha = c(0.05, 0.01),
                 mode="batch",
                 dist.fcts="euclidean")


plot(som.model)
plot(som.model, type="changes")
plot(som.model, type="counts",palette.name=coolBlueHotRed)
plot(som.model, type="dist.neighbours", 
     palette.name=coolBlueHotRed)
plot(som.model, type="mapping", main = "",shape="straight"
     , labels=rownames(dt.scaled), cex=0.45)

#cluster the SOM network with K-Means
codes<- getCodes(som.model)
set.seed(123)
wcss<-rep(0,10)
for (i in 1:10) {
  kmeans.mod<- kmeans(codes, i, nstart=10)
  wcss[i]<-kmeans.mod$tot.withinss
}

plot(1:10, wcss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Within-Cluster Sum of Squares",
     main = "WCSS for Different Numbers of Clusters")

opt.k <- 3
som.clusters <- kmeans(codes, opt.k, nstart=10)
data.clusters <- som.clusters$cluster[som.model$unit.classif]

plot(som.model, type = "mapping", main = "Cluster Mapping",
     bgcol = rainbow(opt.k)[som.clusters$cluster])
add.cluster.boundaries(som.model, som.clusters$cluster)

columns<-list()
for (i in 1:opt.k) {
  #name<-cat("Cluster", i)
  columns[[i]]<-colMeans(dt.dummies[data.clusters == i, ])
}
options(scipen = 999)
summary.tab<- round(data.frame(C1=columns[[1]], 
                               C2=columns[[2]], 
                               C3=columns[[3]]), 2)
summary.tab
xtable(summary.tab)

par(mfrow=c(4,3), cex.main=2, xpd=NA)
for (i in c(1,2,3,4,5,6,10,11,7,8,9,12)) {
  plot(som.model, type = "property", property = codes[, i],
       main = colnames(dt.scaled)[i],
       palette.name=coolBlueHotRed)
}


#######
#Prediction on Response variable
#
set.seed(123)
library(ROCR)
training.index<-sample(1:nrow(dt), 0.7*nrow(dt), replace = F)
train<- dt[training.index,]
test<- dt[- training.index,]

train$Response<- as.factor(train$Response)

index.1<- which(train$Response==1)
index.0<- which(train$Response==0)

#LOGISTIC REGRESSION MODEL
#with oversampling
m<- length(index.0)
boot.ind.1<- sample(index.1, m, replace=T)
boot.ind.0<- sample(index.0, m, replace=F)
train.oversmpl<- train[c(boot.ind.0, boot.ind.1),]

fit<- glm(Response ~ ., data=train.oversmpl[,-c(1)], family="binomial")

pred.probs<- pred.glm<- predict(fit, test[,-c(29)], type="response")
pred<- ifelse(pred.probs>0.5, 1, 0)
acc<- mean(pred==test$Response)
acc
cm.glm<-caret::confusionMatrix(as.factor(pred), as.factor(test$Response),
                       positive="1")

roc.glm <- prediction(predictions = pred.probs, 
                       labels = test$Response)
perf.roc.glm <- performance(roc.glm, measure = "tpr", x.measure = "fpr")
perf.auc.glm <- performance(roc.glm, measure = "auc")
ROC.df.glm <- data.frame(unlist(perf.roc.glm@x.values),
                          unlist(perf.roc.glm@y.values))
colnames(ROC.df.glm) <- c("Specificity","Sensitivity")

pred.probs.glm<- pred.probs #save for distribution plot

#CLASSIFICATION TREE
library(rpart)
cart.model<- rpart(Response ~ ., data=train[,-c(1)],
                   method="class", 
                   parms= list(prior=c(length(index.0)/nrow(train),
                                       length(index.1)/nrow(train)),
                               split="information"))
pred<- predict(cart.model, test[,-c(29)], type="class")
acc<- mean(pred==test$Response)
acc
caret::confusionMatrix(as.factor(pred), as.factor(test$Response),
                       positive="1")

#oversampling
cart.model<- rpart(Response ~ ., data=train.oversmpl[,-c(1)],
                   method="class")
pred.s<- predict(cart.model, test[,-c(29)], type="prob")
pred<- predict(cart.model, test[,-c(29)], type="class")
acc<- mean(pred==test$Response)
acc
cm.cart<-caret::confusionMatrix(as.factor(pred), as.factor(test$Response),
                       positive="1")

roc.cart <- prediction(predictions = pred.s[,"1"], 
                            labels = test$Response)
perf.roc.cart <- performance(roc.cart, measure = "tpr", x.measure = "fpr")
perf.auc.cart <- performance(roc.cart, measure = "auc")
ROC.df.cart <- data.frame(unlist(perf.roc.cart@x.values),
                         unlist(perf.roc.cart@y.values))
colnames(ROC.df.cart) <- c("Specificity","Sensitivity")

#RANDOM FOREST
library(randomForest)

rf.mod<- randomForest(Response ~ ., data=train[,-c(1)],
                      do.trace=T,
                      mtry=6,
                      ntree=1000,
                      importance=T)

pred<- predict(rf.mod, test[,-c(29)], type="class")
caret::confusionMatrix(as.factor(pred), as.factor(test$Response),
                       positive="1")

#oversampling
rf.mod<- randomForest(Response ~ ., data=train.oversmpl[,-c(1)],
                      do.trace=T,
                      mtry=6,
                      ntree=1000,
                      importance=T)

pred<- predict(rf.mod, test[,-c(29)], type="class")
pred.s<- pred.rf<- predict(rf.mod, test[,-c(29)], type="prob")
cm.rf<-caret::confusionMatrix(as.factor(pred), as.factor(test$Response),
                       positive="1")

roc.rf <- prediction(predictions = pred.s[,2], 
                       labels = test$Response)
perf.roc.rf <- performance(roc.rf, measure = "tpr", x.measure = "fpr")
perf.auc.rf <- performance(roc.rf, measure = "auc")
ROC.df.rf <- data.frame(unlist(perf.roc.rf@x.values),
                          unlist(perf.roc.rf@y.values))
colnames(ROC.df.rf) <- c("Specificity","Sensitivity")

pred.probs.rf<- pred.s[,2] #save for distribution plot

#BOOSTING
library(gbm)
ntree.ada = 1000
ada.mod <- gbm(Response ~ ., 
                 data = dt.dummies[training.index,-1], 
                 distribution = "adaboost", 
                 n.trees = ntree.ada,
                 interaction.depth = 10,
                 bag.fraction = 0.5,
                 train.fraction = 0.8,
                 verbose = TRUE,
                 n.cores = 4
)

pred.s <- predict(ada.mod, dt.dummies[-c(training.index),],
                n.trees = ntree.ada, type = 'response')
pred<- ifelse(pred.s>0.5, 1, 0)
cm.ada<- caret::confusionMatrix(as.factor(pred), as.factor(test$Response),
                       positive="1")

roc.ada <- prediction(predictions = pred.s, 
                            labels = test$Response)
perf.roc.ada <- performance(roc.ada, measure = "tpr", x.measure = "fpr")
perf.auc.ada <- performance(roc.ada, measure = "auc")
ROC.df.ada <- data.frame(unlist(perf.roc.ada@x.values),
                         unlist(perf.roc.ada@y.values))
colnames(ROC.df.ada) <- c("Specificity","Sensitivity")

xline <- seq(0,1,0.02)
yline <- seq(0,1,0.02)
xyline <- data.frame(xline,yline)
ggplot() + 
  geom_line(data=ROC.df.glm, aes(x=Specificity, y=Sensitivity, color="GLM")) +
  geom_line(data=ROC.df.cart, aes(x=Specificity, y=Sensitivity, color="CART")) +
  geom_line(data=ROC.df.rf, aes(x=Specificity, y=Sensitivity, color="RF")) +
  geom_line(data=ROC.df.ada, aes(x=Specificity, y=Sensitivity, color="AdaBoost")) +
  geom_line(data=xyline, aes(x=xline, y=yline), color='black',linetype = "dashed") +
  xlab("FPR") + ylab("TPR") +
  scale_colour_manual("Models",
                      values=c("GLM"="blue","CART"="yellow","RF"="red",
                               "AdaBoost"="green")) +
  theme_bw()
#  ggtitle("ROC")

perf.auc.glm@y.values[[1]]
perf.auc.cart@y.values[[1]]
perf.auc.rf@y.values[[1]]
perf.auc.ada@y.values[[1]]

cm.glm
cm.cart
cm.rf
cm.ada

sum.tab.mods<- data.frame(Model=c("Logistic Regression", "CART", "Random Forest", "Boosting"),
                          AUC=c(perf.auc.glm@y.values[[1]],
                                perf.auc.cart@y.values[[1]],
                                perf.auc.rf@y.values[[1]],
                                perf.auc.ada@y.values[[1]]),
                          Acc= c(cm.glm$overall[1], cm.cart$overall[1], cm.rf$overall[1], cm.ada$overall[1]),
                          Sens= c(cm.glm$byClass[1],cm.cart$byClass[1],cm.rf$byClass[1],cm.ada$byClass[1]),
                          Spec= c(cm.glm$byClass[2],cm.cart$byClass[2],cm.rf$byClass[2],cm.ada$byClass[2]))
sum.tab.mods[,2:5]<- round(sum.tab.mods[,2:5], 3)
xtable(sum.tab.mods)

#Distribution plot for predicted prob of response
#rf
p1<-ggplot(data=data.frame(prediction.probability=pred.rf[,2], Response=as.factor(test$Response)),
       aes(x=prediction.probability, group=Response, fill=as.factor(Response))) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_bw() +
  scale_fill_discrete(name = "Response") +
  ylim(0,4.5)  +
  ggtitle("Random Forest model") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))


#glm
p2<-ggplot(data=data.frame(prediction.probability=pred.glm, Response=as.factor(test$Response)),
       aes(x=prediction.probability, group=Response, fill=as.factor(Response))) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_bw() +
  scale_fill_discrete(name = "Response") +
  ylim(0,4.5) +
  ggtitle("Logistic regression model") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

cowplot::plot_grid(p1, p2)

#PROFIT CURVE and PROFIT MAXIMIZATION
#

testnew<-test
cost<- mean(dt$Z_CostContact)
revenue<- mean(dt$Z_Revenue)

profit<- test$Response*(revenue-cost) - (1-test$Response)*cost
testnew$profit<-test$Response*(revenue-cost) - (1-test$Response)*cost
testnew$rf.prediction<-unlist(roc.rf@predictions)
testnew$glm.prediction<-unlist(roc.glm@predictions)

library(dplyr)
library(caret)
testnew<- testnew %>% mutate(rf.rank = rank(desc(rf.prediction),
                                            ties.method = "first"))
testnew<- testnew %>% mutate(glm.rank = rank(desc(glm.prediction),
                                            ties.method = "first"))

testnew <- transform(testnew[order(testnew$glm.rank, decreasing = F), ],
                     cumsum.glm = ave(profit, FUN=cumsum))
testnew <- transform(testnew[order(testnew$rf.rank, decreasing = F), ],
                     cumsum.rf = ave(profit, FUN=cumsum))

testnew$glm.rank.perc<- testnew$glm.rank /max(testnew$glm.rank)
testnew$rf.rank.perc<- testnew$rf.rank /max(testnew$rf.rank)

glm.max.ind<-testnew$glm.rank.perc[which.max(testnew$cumsum.glm)] #0.1704
rf.max.ind<-testnew$rf.rank.perc[which.max(testnew$cumsum.rf)] #0.2272
testnew$cumsum.glm[which(testnew$glm.rank.perc==glm.max.ind)] #417
testnew$cumsum.rf[which(testnew$rf.rank.perc==rf.max.ind)] #446

ggplot() + 
  geom_line(data=testnew, aes(x=glm.rank.perc, y=cumsum.glm, color="GLM"))+
  geom_line(data=testnew, aes(x=rf.rank.perc, y=cumsum.rf, color="RF")) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_vline(xintercept = glm.max.ind, color="blue", linetype="dashed") +
  geom_vline(xintercept = rf.max.ind, color="red", linetype="dashed") +
  xlab("instances") + ylab("profit") +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_manual("Models",
                      values=c("GLM"="blue","RF"="red")) +
  #ggtitle("Cumulative profit") +
  theme_bw()
