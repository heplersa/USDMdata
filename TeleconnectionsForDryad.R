## This script extracts raw teleconnections data from
## https://www.cpc.ncep.noaa.gov
## and formats it to weekly data to match the USDM drought data

## 1. Load libraries
library(stringr)

## 2. ENSO
enso <- read.fwf(
  file=url("https://www.cpc.ncep.noaa.gov/data/indices/wksst9120.for"),
  skip=4,
  widths=c(12, 7, 4, 9, 4, 9, 4, 9, 4)
)

enso <- enso[,c(1,7)]
head(enso)
str_sub(enso[,1], start=-6)

enso[947,] ## earliest date that we need
enso <- enso[-c(1:946),] ## Remove dates before Oct 20, 1999
enso[1185,]
enso <- enso[1:1185,]
head(enso)
tail(enso)

plot(enso[,2], type="l", xlab="Week")

year = as.numeric(str_sub(enso[,1], start=-6, end=-3))
day = as.numeric(str_sub(enso[,1], start=2, end=3))

mon.txt = str_sub(enso[,1], start=4, end=6)
mon = array()
mon[mon.txt=="JAN"] = "01"
mon[mon.txt=="FEB"] = "02"
mon[mon.txt=="MAR"] = "03"
mon[mon.txt=="APR"] = "04"
mon[mon.txt=="MAY"] = "05"
mon[mon.txt=="JUN"] = "06"
mon[mon.txt=="JUL"] = "07"
mon[mon.txt=="AUG"] = "08"
mon[mon.txt=="SEP"] = "09"
mon[mon.txt=="OCT"] = "10"
mon[mon.txt=="NOV"] = "11"
mon[mon.txt=="DEC"] = "12"

enso.wk.pre = data.frame(year, mon, day, enso[,2])
head(enso.wk.pre)
tail(enso.wk.pre)
## Matches by one day (i.e. 1/5/2000 goes with 1/4/2000 in final dataset)

dim(enso.wk.pre)
enso.wk = enso.wk.pre[12:1185,]
head(enso.wk)
enso.14 = enso.wk[,4]
enso.28 = enso.wk[,4]
enso.84 = enso.wk[,4]
for (t in 1:length(enso.14)){
  enso.14[t] = mean(enso.wk.pre[((t+11)-1):(t+11),4])
  enso.28[t] = mean(enso.wk.pre[((t+11)-3):(t+11),4])
  enso.84[t] = mean(enso.wk.pre[((t+11)-11):(t+11),4])
}
enso.wk = cbind(enso.wk, enso.14, enso.28, enso.84)
head(enso.wk)

## 3. PNA
pna <- read.fwf(
  file=url("https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.pna.index.b500101.current.ascii"),
  skip=0,
  widths=c(4, 3, 3, 7)
)
colnames(pna) <- c("year","month","day","pna")

wh = which(pna$year == 2000 & pna$month == 1 & pna$day == 4)
## upper=ceiling((dim(pna)[1] - wh)/7)
upper = 1174

pna.wk = pna[(wh+7*(seq(1,upper)-1)),]
tail(pna.wk)

pna.avg = rowMeans(matrix(as.numeric(pna$pna[(wh-6):26477]), ncol = 7, byrow=TRUE), na.rm=TRUE)
tail(pna.avg)
pna.wk[,4] <- pna.avg


## 4. NAO
nao <- read.fwf(
  file=url("https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.nao.index.b500101.current.ascii"),
  skip=0,
  widths=c(4, 3, 3, 7)
)
colnames(nao) <- c("year","month","day","nao")

wh = which(nao$year == 2000 & nao$month == 1 & nao$day == 4)
nao.wk = nao[(wh+7*(seq(1,upper)-1)),]
tail(nao.wk)

nao.avg = rowMeans(matrix(as.numeric(nao$nao[(wh-6):26477]), ncol = 7, byrow=TRUE), na.rm=TRUE)
tail(nao.avg)
nao.wk[,4] <- nao.avg


## 5. AO
ao <- read.fwf(
  file=url("https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.ao.index.b500101.current.ascii"),
  skip=0,
  widths=c(4, 3, 3, 7)
)
colnames(ao) <- c("year","month","day","ao")

wh = which(ao$year == 2000 & ao$month == 1 & ao$day == 4)
ao.wk = ao[(wh+7*(seq(1,upper)-1)),]
tail(ao.wk)

ao.avg = rowMeans(matrix(as.numeric(ao$ao[(wh-6):26477]), ncol = 7, byrow=TRUE), na.rm=TRUE)
tail(ao.avg)
ao.wk[,4] <- ao.avg

## 6. Combine teleconnections into common format
head(pna.wk)
time = pna.wk$day*1 + pna.wk$month*100 + pna.wk$year*10000
tele <- data.frame(time, pna.wk$pna, nao.wk$nao, ao.wk$ao)
enso.wk[1:dim(tele)[1],4:7]
tele <- cbind(tele, enso.wk[1:dim(tele)[1],4:7])
colnames(tele) <- c("time","pna","nao","ao","enso","enso14","enso28","enso84")

tele$pna <- round(tele$pna, 2)
tele$nao <- round(tele$nao, 2)
tele$ao <- round(tele$ao, 2)
tele$enso <- round(tele$enso, 2)
tele$enso14 <- round(tele$enso14, 3)
tele$enso28 <- round(tele$enso28, 3)
tele$enso84 <- round(tele$enso84, 3)

head(tele)
write.csv(file="Teleconnections/tele3.csv", tele, row.names=FALSE)

