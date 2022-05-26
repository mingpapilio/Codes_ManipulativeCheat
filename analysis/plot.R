library(data.table)
##
id_start<- 1
id_end<- 1
##
raw<-fread("zlog.txt", sep="\t")
## Data Processing
colnames(raw)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10","NA")
raw<-raw[,!"NA"]
leng<- length(raw[,T])
## Plotting
dev.new()
par(mar=c(5, 5, 2, 2))
plot(c(0,0), xlim=c(0, 3000), ylim=c(0,1), col="white", xlab="Generation", ylab="Trait values")
# plot(c(0,0), xlim=c(0, raw[leng,T]), ylim=c(0,1), col="white", xlab="Generation", ylab="Trait values")
for (i in id_start: id_end){
  idx<- paste("rep_", i, sep="")
  points(raw[1:600,T], raw[1:600,get(idx)], type="l", col=rgb(228/255,0,127/255,1),lwd=1)
}
##
raw<-fread("x1log.txt", sep="\t")
## Data Processing
colnames(raw)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10","NA")
raw<-raw[,!"NA"]
for (i in id_start: id_end){
  idx<- paste("rep_", i, sep="")
  points(raw[1:600,T], raw[1:600,get(idx)], type="l", col=rgb(46/255,167/255,224/255,1),lwd=1)
}
##
raw<-fread("u1log.txt", sep="\t")
## Data Processing
colnames(raw)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10","NA")
raw<-raw[,!"NA"]
for (i in id_start: id_end){
  idx<- paste("rep_", i, sep="")
  points(raw[1:600,T], raw[1:600,get(idx)], type="l", col=rgb(248/255,182/255,45/255,1),lwd=1)
}
##
raw<-fread("x2log.txt", sep="\t")
## Data Processing
colnames(raw)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10","NA")
raw<-raw[,!"NA"]
for (i in id_start: id_end){
  idx<- paste("rep_", i, sep="")
  points(raw[1:600,T], raw[1:600,get(idx)], type="l", col=rgb(23/255,42/255,136/255,1),lwd=1)
}
##
raw<-fread("u2log.txt", sep="\t")
## Data Processing
colnames(raw)<- c("T", "rep_1", "rep_2", "rep_3", "rep_4", "rep_5", "rep_6", "rep_7", "rep_8", "rep_9", "rep_10","NA")
raw<-raw[,!"NA"]
for (i in id_start: id_end){
  idx<- paste("rep_", i, sep="")
  points(raw[1:600,T], raw[1:600,get(idx)], type="l", col=rgb(234/255,85/255,20/255,1),lwd=1)
}
##
