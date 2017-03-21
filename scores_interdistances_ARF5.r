(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)
source("PNfonctions.r")                 # fonctions auxiliaires

#-------------------------------------read matrices ------------------------------------------


pfm_ARF5<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF5 <- round((t(as.matrix(pfm_ARF5)))*nRegion)+1 ;pfm_ARF5
pfm_ARF5 <- cbind(rep(151,4),rep(151,4),pfm_ARF5)
maxi_ARF5 <- apply(pfm_ARF5,FUN=max, 2)
maxi_ARF5 <- matrix(nrow=4, rep(maxi_ARF5,4),byrow=TRUE)
pwm_ARF5 <- log(pfm_ARF5/maxi_ARF5)
pwm_ARF5_rev <- pwm_ARF5 - minScore(pwm_ARF5)/dim(pwm_ARF5)[2] 
pwm_ARF5 <-  reverseComplement(pwm_ARF5_rev) ; pwm_ARF5

#-------------------------------------read fasta files-----------------------------------------


ARF5_pos <- readDNAStringSet('MP_pos.fas')
ARF5_neg <- readDNAStringSet('ARF2_pos.fas')#_wo_ARF2.fas') 
width_pos <- width(ARF5_pos)
width_neg <- width(ARF5_neg)

seq_pos <- as.character(ARF5_pos)
seq_neg <- as.character(ARF5_neg)


#-------------------------------------Compute Scores-----------------------------------------

th <- maxScore(pwm_ARF5) - 12

#pos

scores_ARF5_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF5)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF5))

scores_ARF5_rev_pos <- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF5_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF5_rev))

density_pos <- (sapply(FUN=sum,lapply(FUN=">",scores_ARF5_pos,th))+sapply(FUN=sum,lapply(FUN=">",scores_ARF5_rev_pos,th)))/(width_pos*2)

#neg

scores_ARF5_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF5)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF5))

scores_ARF5_rev_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF5_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF5_rev))

density_neg <- (sapply(FUN=sum,lapply(FUN=">",scores_ARF5_neg,th))+sapply(FUN=sum,lapply(FUN=">",scores_ARF5_rev_neg,th)))/(width_neg*2)

#-----------------------------------Compute Scores Interdistances-----------------------------

#initialise

scores_DR_pos<- matrix(0,nrow=length(ARF5_pos),ncol=21)
scores_DR_pos_rev<- matrix(0,nrow=length(ARF5_pos),ncol=21)
scores_ER_pos<- matrix(0,nrow=length(ARF5_pos),ncol=21)
scores_IR_pos<- matrix(0,nrow=length(ARF5_pos),ncol=21)

scores_DR_neg<- matrix(0,nrow=length(ARF5_neg),ncol=21)
scores_DR_neg_rev<- matrix(0,nrow=length(ARF5_neg),ncol=21)
scores_ER_neg<- matrix(0,nrow=length(ARF5_neg),ncol=21)
scores_IR_neg<- matrix(0,nrow=length(ARF5_neg),ncol=21)


# Compute DR 

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_DR_pos[i,j-6] <- max(scores_ARF5_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5)[2])] + scores_ARF5_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5)[2])])
    }
}


for(i in 1:length(ARF5_neg))
{
    for (j in 7:27)
    {
        scores_DR_neg[i,j-6] <- max(scores_ARF5_neg[[i]][j:(width_neg[i]-dim(pwm_ARF5)[2])] + scores_ARF5_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF5)[2])])
    }
}

# Compute DR rev

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_DR_pos_rev[i,j-6] <- max(scores_ARF5_rev_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5_rev)[2])] + scores_ARF5_rev_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5_rev)[2])])
    }
}


for(i in 1:length(ARF5_neg))
{
    for (j in 7:27)
    {
        scores_DR_neg_rev[i,j-6] <- max(scores_ARF5_rev_neg[[i]][j:(width_neg[i]-dim(pwm_ARF5_rev)[2])] + scores_ARF5_rev_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF5_rev)[2])])
    }
}

# Best score between DR and DR rev

scores_DR_pos <- ifelse(scores_DR_pos > scores_DR_pos_rev, scores_DR_pos,scores_DR_pos_rev)
scores_DR_neg <- ifelse(scores_DR_neg > scores_DR_neg_rev, scores_DR_neg,scores_DR_neg_rev)

# Compute IR 

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_IR_pos[i,j-6] <- max(scores_ARF5_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5)[2])] + scores_ARF5_rev_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5_rev)[2])])
    }
}


for(i in 1:length(ARF5_neg))
{
    for (j in 7:27)
    {
        scores_IR_neg[i,j-6] <- max(scores_ARF5_neg[[i]][j:(width_neg[i]-dim(pwm_ARF5)[2])] + scores_ARF5_rev_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF5_rev)[2])])
    }
}

# Compute ER

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_ER_pos[i,j-6] <- max(scores_ARF5_rev_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5_rev)[2])] + scores_ARF5_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5)[2])])
    }
}


for(i in 1:length(ARF5_neg))
{
    for (j in 7:27)
    {
        scores_ER_neg[i,j-6] <- max(scores_ARF5_rev_neg[[i]][j:(width_neg[i]-dim(pwm_ARF5_rev)[2])] + scores_ARF5_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF5)[2])])
    }
}



#------------------------Meilleurs scores--------------------------------
pos <- apply(FUN=max, cbind(scores_DR_pos,scores_IR_pos,scores_ER_pos),MARGIN=1)
neg <- apply(FUN=max, cbind(scores_DR_neg,scores_IR_neg,scores_ER_neg),MARGIN=1)


scores_ER_pos2 <- scores_ER_pos
scores_ER_neg2 <- scores_ER_neg
scores_IR_pos2 <- scores_IR_pos
scores_IR_neg2 <- scores_IR_neg
scores_DR_pos2 <- scores_DR_pos
scores_DR_neg2 <- scores_DR_neg
scores_DR_pos2[,5] <- scores_DR_pos2[,5] + 3
scores_DR_neg2[,5] <- scores_DR_neg2[,5] + 3
scores_DR_pos2[,6] <- scores_DR_pos2[,6] + 2
scores_DR_neg2[,6] <- scores_DR_neg2[,6] + 2
scores_IR_pos2[,14] <- scores_IR_pos2[,14] + 6
scores_IR_neg2[,14] <- scores_IR_neg2[,14] + 6
scores_IR_pos2[,1] <- scores_IR_pos2[,1] + 9
scores_IR_neg2[,1] <- scores_IR_neg2[,1] + 9
scores_IR_pos2[,3] <- scores_IR_pos2[,3] + 6
scores_IR_neg2[,3] <- scores_IR_neg2[,3] + 6
scores_IR_pos2[,13] <- scores_IR_pos2[,13] + 2
scores_IR_neg2[,13] <- scores_IR_neg2[,13] + 2
pos_pen <- apply(FUN=max, cbind(scores_DR_pos2,scores_IR_pos2,scores_ER_pos2),MARGIN=1)
neg_pen <- apply(FUN=max, cbind(scores_DR_neg2,scores_IR_neg2,scores_ER_neg2),MARGIN=1)
pos_pen2 <- pos_pen + density_pos *100
neg_pen2 <- neg_pen + density_neg *100
rc1 = ROCcurve(pos,neg) # fait la roc
X <- rc1$XY[1,]
Y <- rc1$XY[2,]
AU <- rc1$AUC
A <- as.character(round(AU,4))
AUC <- paste("AUC = ", A,sep="")
rc2 = ROCcurve(pos_pen,neg_pen) # fait la roc
X_pen <- rc2$XY[1,]
Y_pen <- rc2$XY[2,]
AU_pen<- rc2$AUC
A_pen<- as.character(round(AU_pen,4))
AUC_pen<- paste("AUC with model = ", A_pen,sep="")
rc3 = ROCcurve(pos_pen2,neg_pen2) # fait la roc
X_pen2 <- rc3$XY[1,]
Y_pen2 <- rc3$XY[2,]
AU_pen2<- rc3$AUC
A_pen2<- as.character(round(AU_pen2,4))
AUC_pen2<- paste("AUC with penalties and density = ", A_pen2,sep="")
{plot(X,Y,type="l",col="red",lwd=2,
      ylab="ARF5",xlab="ARF2",
      main="ARF5 vs ARF2")}
lines(X_pen,Y_pen,col='cornflowerblue',lwd=2)
lines(X_pen2,Y_pen2,col='green4',lwd=2)
{legend(0.3,0.15,legend=c(AUC,AUC_pen,AUC_pen2),
        col=c("red","cornflowerblue","green4"),lty=rep(1,4),lwd=rep(2,4))}
dev.copy(device = png, filename = 'ROC_ARF5_vs_ARF2.png', width = 800, height = 600)
dev.off()



## pn01 <- c(rep(0,length(ARF2_pos)),rep(1,length(ARF2_neg)))
    
## lm1 <- glm(pn01~rbind(scores_ER_pos,scores_ER_neg),family=binomial)   # variable 1

## an1 <- anova(lm1,test="Chisq")
