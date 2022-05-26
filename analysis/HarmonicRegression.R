library(data.table)
library(HarmonicRegression)
## Read the file
raw<-fread("zlog.txt", sep="\t")
## Data Processing if it detects 12 columns
colnames(raw)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10","NA")
raw<-raw[,!"NA"]
## Harmonic regression
raw<- as.matrix(raw)
leng<- 360
pars<- matrix(NA, leng, 11)
cis<- matrix(NA, leng,11)
colnames(pars)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10")
colnames(cis)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10")
for(i in 12:360){
  pars[i,1]<- 5*i
  cis[i,1]<- 5*i
  pars[i,2:11]<- harmonic.regression(inputts=raw[150:600,2:11], inputtime=raw[150:600, 1], Tau= i)$pars$amp
  cis[i,2:11]<- harmonic.regression(inputts=raw[150:600,2:11], inputtime=raw[150:600, 1], Tau= i)$ci$amp
}
## Taking averages over the repeats
log<- matrix(NA, leng, 3)
colnames(log)<- c("T", "mean","ci")
for(i in 1:leng){
  log[i,1]<- pars[i,1]
  log[i,2]<- mean(pars[i,2:11])
  log[i,3]<- mean(cis[i,2:11])
}
## Plotting
dev.new()
par(mar=c(5, 5, 2, 2))
plot(c(1,1),col="white",xlim=c(60,1800),ylim=c(0,0.05),xlab="Period (generation)", ylab="Amplitude")
for(i in 1:leng){
  segments(x0=log[i,1],y0=(log[i,2]-log[i,3]),y1=(log[i,2]+log[i,3]),col="gray",lwd=2)
}
points(log[,1], log[,2], type="l")
