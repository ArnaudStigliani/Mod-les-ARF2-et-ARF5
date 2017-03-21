rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)
source("PNfonctions.r")                 # fonctions auxiliaires

#-------------------------------------read matrices ------------------------------------------


pfm_ARF2<- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF2 <- round((t(as.matrix(pfm_ARF2)))*nRegion)+1 ;pfm_ARF2
maxi_ARF2 <- apply(pfm_ARF2,FUN=max, 2)
maxi_ARF2 <- matrix(nrow=4, rep(maxi_ARF2,4),byrow=TRUE)
pwm_ARF2 <- log(pfm_ARF2/maxi_ARF2)
pwm_ARF2_rev <- pwm_ARF2 - minScore(pwm_ARF2)/dim(pwm_ARF2)[2] 

pwm_ARF2 <-  reverseComplement(pwm_ARF2_rev) ; pwm_ARF2

#-------------------------------------read fasta files-----------------------------------------


ARF2_pos <- readDNAStringSet('ARF2_pos.fas')[(1:1000)]
ARF2_neg <- readDNAStringSet('ARF2_negativesetbuilder.fas')[(1:1000)]
width_pos <- width(ARF2_pos)
width_neg <- width(ARF2_neg)

seq_pos <- as.character(ARF2_pos)
seq_neg <- as.character(ARF2_neg)


#-------------------------------------Compute Scores-----------------------------------------

#pos

scores_ARF2_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF2)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2))

scores_ARF2_rev_pos <- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF2_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2_rev))

#neg

scores_ARF2_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF2)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2))

scores_ARF2_rev_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF2_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2_rev))

#-----------------------------------Compute Scores Interdistances-----------------------------

#initialise

scores_DR_pos<- matrix(0,nrow=length(ARF2_pos),ncol=21)
scores_DR_pos_rev<- matrix(0,nrow=length(ARF2_pos),ncol=21)
scores_ER_pos<- matrix(0,nrow=length(ARF2_pos),ncol=21)
scores_IR_pos<- matrix(0,nrow=length(ARF2_pos),ncol=21)

scores_DR_neg<- matrix(0,nrow=length(ARF2_neg),ncol=21)
scores_DR_neg_rev<- matrix(0,nrow=length(ARF2_neg),ncol=21)
scores_ER_neg<- matrix(0,nrow=length(ARF2_neg),ncol=21)
scores_IR_neg<- matrix(0,nrow=length(ARF2_neg),ncol=21)


# Compute DR 

for(i in 1:length(ARF2_pos))
{
    for (j in 7:27)
    {
        scores_DR_pos[i,j-6] <- max(scores_ARF2_pos[[i]][j:(width_pos[i]-dim(pwm_ARF2)[2])] + scores_ARF2_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF2)[2])])
    }
}


for(i in 1:length(ARF2_neg))
{
    for (j in 7:27)
    {
        scores_DR_neg[i,j-6] <- max(scores_ARF2_neg[[i]][j:(width_neg[i]-dim(pwm_ARF2)[2])] + scores_ARF2_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF2)[2])])
    }
}

# Compute DR rev

for(i in 1:length(ARF2_pos))
{
    for (j in 7:27)
    {
        scores_DR_pos_rev[i,j-6] <- max(scores_ARF2_rev_pos[[i]][j:(width_pos[i]-dim(pwm_ARF2_rev)[2])] + scores_ARF2_rev_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF2_rev)[2])])
    }
}


for(i in 1:length(ARF2_neg))
{
    for (j in 7:27)
    {
        scores_DR_neg_rev[i,j-6] <- max(scores_ARF2_rev_neg[[i]][j:(width_neg[i]-dim(pwm_ARF2_rev)[2])] + scores_ARF2_rev_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF2_rev)[2])])
    }
}

# Best score between DR and DR rev

scores_DR_pos <- ifelse(scores_DR_pos > scores_DR_pos_rev, scores_DR_pos,scores_DR_pos_rev)
scores_DR_neg <- ifelse(scores_DR_neg > scores_DR_neg_rev, scores_DR_neg,scores_DR_neg_rev)

# Compute IR 

for(i in 1:length(ARF2_pos))
{
    for (j in 7:27)
    {
        scores_IR_pos[i,j-6] <- max(scores_ARF2_pos[[i]][j:(width_pos[i]-dim(pwm_ARF2)[2])] + scores_ARF2_rev_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF2_rev)[2])])
    }
}


for(i in 1:length(ARF2_neg))
{
    for (j in 7:27)
    {
        scores_IR_neg[i,j-6] <- max(scores_ARF2_neg[[i]][j:(width_neg[i]-dim(pwm_ARF2)[2])] + scores_ARF2_rev_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF2_rev)[2])])
    }
}

# Compute ER

for(i in 1:length(ARF2_pos))
{
    for (j in 7:27)
    {
        scores_ER_pos[i,j-6] <- max(scores_ARF2_rev_pos[[i]][j:(width_pos[i]-dim(pwm_ARF2_rev)[2])] + scores_ARF2_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF2)[2])])
    }
}


for(i in 1:length(ARF2_neg))
{
    for (j in 7:27)
    {
        scores_ER_neg[i,j-6] <- max(scores_ARF2_rev_neg[[i]][j:(width_neg[i]-dim(pwm_ARF2_rev)[2])] + scores_ARF2_neg[[i]][1:(width_neg[i]-j+1-dim(pwm_ARF2)[2])])
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
scores_ER_pos2[,8] <- scores_ER_pos2[,8] + 7
scores_ER_neg2[,8] <- scores_ER_neg2[,8] + 7
scores_ER_pos2[,9] <- scores_ER_pos2[,9] + 7
scores_ER_neg2[,9] <- scores_ER_neg2[,9] + 7
## scores_IR_pos2[,12] <- scores_IR_pos2[,12] +3
## scores_IR_neg2[,12] <- scores_IR_neg2[,12] +3
## scores_IR_pos2[,10] <- scores_IR_pos2[,10] +3
## scores_IR_neg2[,10] <- scores_IR_neg2[,10] +3
scores_DR_pos2 <- scores_DR_pos2 - 6
scores_DR_neg2 <- scores_DR_neg2 - 6
pos_pen <- apply(FUN=max, cbind(scores_DR_pos2,scores_IR_pos2,scores_ER_pos2),MARGIN=1)
neg_pen <- apply(FUN=max, cbind(scores_DR_neg2,scores_IR_neg2,scores_ER_neg2),MARGIN=1)
rc1 = ROCcurve(pos,neg) # fait la roc
X <- rc1$XY[1,]
Y <- rc1$XY[2,]
rc2 = ROCcurve(pos_pen,neg_pen) # fait la roc
X_pen <- rc2$XY[1,]
Y_pen <- rc2$XY[2,]
AU <- rc1$AUC
A <- as.character(round(AU,4))
AUC <- paste("AUC = ", A,sep="")
AU_pen<- rc2$AUC
A_pen<- as.character(round(AU_pen,4))
AUC_pen<- paste("AUC with penalties = ", A_pen,sep="")
{plot(X,Y,type="l",col="red",lwd=2,
      ylab="ARF2",xlab="negative set",
      main="ARF2 vs negative set")}
lines(X_pen,Y_pen,col='cornflowerblue',lwd=2)
{legend(0.5,0.3,legend=c(AUC,AUC_pen),
        col=c("red","cornflowerblue"),lty=rep(1,4),lwd=rep(2,4))}
dev.copy(device = png, filename = 'ROC_with_penalties_1000.png', width = 800, height = 600) 
dev.off()



## pn01 <- c(rep(0,length(ARF2_pos)),rep(1,length(ARF2_neg)))
    
## lm1 <- glm(pn01~rbind(scores_ER_pos,scores_ER_neg),family=binomial)   # variable 1

## an1 <- anova(lm1,test="Chisq")
