################################################
#linear model analysis                         #
#values: "height and diameter"                 #
#Fixed factors: "genotype" & "time"            #
#                                              #
################################################

#load the library "lme4"#
library(lme4)

#data input#
input.file <- file.choose()
dat <- read.delim(input.file, header=T)

#define fixed factors
genotype <- levels(dat$genotype)
time <- levels(dat$time)

#Normarity test by Q-Q plot
par(mfrow=c(1,2))
qqnorm(dat$heigh)
qqline(dat$height, col="red")

#mixed linear model
dat_lmer <- summary( lmer( height ~ genotype:time - 1 + (1|rep), data=dat))

#Linearity, Homoscedasticity and Normatilty test
par(mfrow=c(2,2))
plot(dat[,"height"], dat_lmer$residuals)
abline(0,0,col="red")
qqnorm(dat_lmer$residuals)
qqline(dat_lmer$residuals, col="red")

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_diameter_week1_5.txt", sep="\t", row.names=T,col.names=T, quote=F)
#SE of fixed effects = sqrt(diag(vcov(dat_lmer)))

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtime interaction
genotype.time <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.timecomp <- c()								
for (i in 1:(length(genotype.time)-1)){							
  for (j in (i+1):length(genotype.time)){
    genotype.timecomp <- c(genotype.timecomp, paste(genotype.time[i], genotype.time[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.timecomp)
genotype.treatmentcomp = gsub("time", "", genotype.timecomp)


#Calculate p-values (genotype:treatment vs genotype:time)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.time) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.timecomp
write.table(p.val, "p.val_all_diameter_week1_5.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.01, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_diameter.txt", sep="\t", row.names=T,col.names=F, quote=F)

 #Calculate p-values (Diff of Diff)

geno.no <- length(genotype)
time.no <- length(time)
a.genotype <- c(1:geno.no)

Effect <- estim[a.genotype] - estim[a.genotype + geno.no]
names(Effect) <- genotype

id.mat <- matrix( 0, ncol=1, nrow = ncol( vcov ) )  # all 0, 2*geno.num x 1 matrix
p.val <- c()
p.names <- c()
for (i in 1:(geno.no - 1) ) {
  for ( j in (i + 1):geno.no ) {
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- id.mat.x[ j + geno.no, 1 ] <- 1
    id.mat.x[ i + geno.no, 1 ] <- id.mat.x[ j, 1 ] <- -1
    se <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    dif <- Effect[i] - Effect[j]
    t.val <- abs(dif) / se
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
    p.names <- c( p.names, paste(genotype[i],genotype[j], sep='-') )
  }
}
names(p.val) <- p.names
write.table(p.val, "p.val_DofD_diameter_week1-2.txt", sep="\t", row.names=T,col.names=F, quote=F)

p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.05, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter.txt", sep="\t", row.names=T,col.names=F, quote=F)



################################################
#linear model analysis                         #
#values: "SCW thickness                        #
#Fixed factors: "genotype"                     #
#                                              #
################################################
library(lme4)

#data input#
input.file <- file.choose()
dat <- read.delim(input.file, header=T)


qqnorm(dat$Area.ZW.inches2)
qqline(dat$Area.ZW.inches2, col='red')

shapiro.test(dat$Area.ZW.inches2)
shapiro.test(rnorm(100))

genotype <- levels(dat$line)
replicate <- levels(dat$replicate)

#Linear model analysis
dat_lmer <- summary(lmer( Area.ZW.inches2 ~ line - 1 + (1|replicate), data=dat))

#Linearity, Homoscedasticity and Normatilty test
par(mfrow=c(2,2))
plot(dat[,"Area.ZW.inches2"], dat_lmer$residuals)
abline(0,0,col="red")
qqnorm(dat_lmer$residuals)
qqline(dat_lmer$residuals, col="red")

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_SCWERF139.txt", sep="\t", row.names=T,col.names=T, quote=F)

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp
write.table(p.val, "p.val_all_SCWERF139.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.001, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_SCWthickness_ERF139.txt", sep="\t", row.names=T,col.names=F, quote=F)


################################################
#linear model analysis                         #
#values: TW:OW ratios                          #
#Fixed factors: "genotype"                     #
#                                              #
################################################
library(lme4)

#data input#
input.file <- file.choose()
dat <- read.delim(input.file, header=T)


qqnorm(dat$TW_OW)
qqline(dat$TW_OW, col='red')

shapiro.test(dat$TW_OW)
shapiro.test(rnorm(100))

genotype <- levels(dat$line)
replicate <- levels(dat$replicate)

#Linear model analysis
dat_lmer <- summary(lmer( TW_OW ~ line - 1 + (1|replicate), data=dat))

#Linearity, Homoscedasticity and Normatilty test
par(mfrow=c(2,2))
plot(dat[,"TW_OW"], dat_lmer$residuals)
abline(0,0,col="red")
qqnorm(dat_lmer$residuals)
qqline(dat_lmer$residuals, col="red")

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_TW_OWERF139SRDX.txt", sep="\t", row.names=T,col.names=T, quote=F)

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp
write.table(p.val, "p.val_all_TW_OW_ERF139SRDX.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.05, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_TWOW_ERF139OE.txt", sep="\t", row.names=T,col.names=F, quote=F)



################################################
#linear model analysis                         #
#values: "Vessel No and size in NW"            #
#Fixed factors: "genotype"                     #
#                                              #
################################################

library(lme4)

#data input#
input.file <- file.choose()
dat <- read.delim(input.file, header=T)


qqnorm(dat$VesselNo_mm2)
qqline(dat$VesselNo_mm2, col='red')

shapiro.test(dat$VesselNo_mm2)
shapiro.test(rnorm(100))

genotype <- levels(dat$genotype)
replicate <- levels(dat$replicate)

#Linear model analysis
dat_lmer <- summary(lmer( VesselNo_mm2 ~ genotype - 1 + (1|replicate), data=dat))

#Linearity, Homoscedasticity and Normatilty test
par(mfrow=c(2,2))
plot(dat[,"VesselNo_mm2"], dat_lmer$residuals)
abline(0,0,col="red")
qqnorm(dat_lmer$residuals)
qqline(dat_lmer$residuals, col="red")

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_VesselNo_NWERF139.txt", sep="\t", row.names=T,col.names=T, quote=F)

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp
write.table(p.val, "p.val_VesselNo_NWERF139.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.01, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_VesselNo_NWERF139.txt", sep="\t", row.names=T,col.names=F, quote=F)



qqnorm(dat$area_um2)
qqline(dat$area_um2, col='red')

shapiro.test(dat$area_um2)
shapiro.test(rnorm(100))

genotype <- levels(dat$genotype)
replicate <- levels(dat$replicate)

#Linear model analysis
dat_lmer <- summary(lmer( area_um2 ~ genotype - 1 + (1|replicate), data=dat))

#Linearity, Homoscedasticity and Normatilty test
par(mfrow=c(2,2))
plot(dat[,"area_um2"], dat_lmer$residuals)
abline(0,0,col="red")
qqnorm(dat_lmer$residuals)
qqline(dat_lmer$residuals, col="red")

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_area_um2_NWERF139.txt", sep="\t", row.names=T,col.names=T, quote=F)

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp
write.table(p.val, "p.val_area_um2_NWERF139.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.01, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_area_um2_NWERF139.txt", sep="\t", row.names=T,col.names=F, quote=F)


################################################
#linear model analysis                         #
#values: "Vessel No and size in OW and TW"     #
#Fixed factors: "genotype" & "treatment"       #
#                                              #
################################################

#load the library "lme4"#
library(lme4)

#data input#
input.file <- file.choose()
dat <- read.delim(input.file, header=T)

#define fixed factors
genotype <- levels(dat$genotype)
treatment <- levels(dat$treatment)

#Normarity test by Q-Q plot
par(mfrow=c(1,2))
qqnorm(dat$area_um2)
qqline(dat$area_um2, col="red")

#mixed linear model
dat_lmer <- summary( lmer( vesselNo ~ genotype:treatment - 1 + (1|replicate), data=dat))

#Linearity, Homoscedasticity and Normatilty test
par(mfrow=c(2,2))
plot(dat[,"vesselNo"], dat_lmer$residuals)
abline(0,0,col="red")
qqnorm(dat_lmer$residuals)
qqline(dat_lmer$residuals, col="red")

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_numbervessels_TWOW_ERF139SRDX.txt", sep="\t", row.names=T,col.names=T, quote=F)
#SE of fixed effects = sqrt(diag(vcov(dat_lmer)))

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)
genotype.treatmentcomp = gsub("treatment", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp
write.table(p.val, "p.val_all_number_vessels_TWOW_ERF139SRDX.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.05, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_area_um2_TWOW_ERF139SRDX.txt", sep="\t", row.names=T,col.names=F, quote=F)


################################################
#linear model analysis                         #
#values: "chemistry in NW"                     #
#Fixed factors: "genotype"                     #
#                                              #
################################################

library(lme4)

#data input#
input.file <- file.choose()
dat <- read.delim(input.file, header=T)


genotype <- levels(dat$genotype)
replicate <- levels(dat$replicate)

#Linear model analysis
dat_lmer <- summary(lmer( L ~ genotype - 1 + (1|replicate), data=dat))

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]
write.table(dat_lmer$coefficients, "estim_L_Pyrolisis.txt", sep="\t", row.names=T,col.names=T, quote=F)

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp
write.table(p.val, "p.val_L_pyrolisis.txt", sep="\t", row.names=T,col.names=F, quote=F)

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.01, Letters=c(letters, LETTERS, "."), reverse=F)
write.table(p.val.lett$Letters, "p.val_letter_L.txt", sep="\t", row.names=T,col.names=F, quote=F)

