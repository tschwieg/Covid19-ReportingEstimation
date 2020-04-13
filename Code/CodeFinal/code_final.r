rm(list = ls())
## We use dplyr for filtering data and aggregating
library( dplyr)
library( ggplot2 )
library( ggrepel)
library( stargazer )

## This gives a time series of infections around the world merged with population data.
worldInfection <- read.csv("../../Data/mergedData.csv")


## Table giving entry into Iceland in 2020 by country of origin
iceLandEntry <- read.csv("../../Data/IcelandEntry.csv" )

## This function returns the dataframe for available states/counties given by cutoffDate
## where alpha can be computed by lm( Confirmed ~ Reg - 1, data = MergedData )
EstimateAlpha <- function( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate,
                          iceLandEntry, excludeEarly, useItaly, useDECODE, Lag, fullEU ){
    stateInfection <- read.csv("../../Data/time_series_covid19_confirmed_US.csv" )

    if( excludeEarly ){

        ## We only exclude Seattle. Consider also the removal of San
        ## Diego and Maricopa AZ as they also had early infection.
        stateInfection <- stateInfection %>% filter( !(Combined_Key %in% c( "King, Washington, US" ) ) )
    }
    stateInfection <- merge( x=statePop[,c("Pop", "Combined_Key")], y = stateInfection,
                            by="Combined_Key" )

    chinaData <- worldInfection[match( "China", worldInfection$Country.Region ),]
    italyData <- worldInfection[match( "Italy", worldInfection$Country.Region ),]
    spainData <- worldInfection[match( "Spain", worldInfection$Country.Region ),]
    germanData <- worldInfection[match( "Germany", worldInfection$Country.Region ),]
    ukData <- worldInfection[match( "United Kingdom", worldInfection$Country.Region ),]

   
    ## Get the dates from the string cutoffDate
    PeriodT1Col <- 0
    if( substr( cutoffDate, 1, 2 ) == "03" ){
        PeriodT1Col <- match( "X3.1.20", names(worldInfection) ) - 1 +
            as.numeric( substr(cutoffDate, 4,5 ) ) 
        T1ColUS <- match( "X3.1.2020", names(stateInfection) ) - 1 +
            as.numeric( substr(cutoffDate, 4,5 ) )
    }else if( substr( cutoffDate, 1, 2 ) == "02" ){
        PeriodT1Col <- match( "X2.1.20", names(worldInfection) ) - 1 +
            as.numeric( substr(cutoffDate, 4,5 ) )
        T1ColUS <- match( "X2.1.2020", names(stateInfection) ) - 1 +
            as.numeric( substr(cutoffDate, 4,5 ) )
    }
    

    ## Note this only works for Feb T0. Need to modify for other months.
    PeriodT0ColChina <- match( "X2.1.20", names(worldInfection) ) - 1 +
        as.numeric( substr(t0DateChina, 4,5 ) )

    ## Italy begins to record infections much later than China, need a separate date
    PeriodT0ColItaly <- match( "X2.1.20", names(worldInfection) ) - 1 +
        as.numeric( substr(t0DateItaly, 4,5 ) )

    
    ## We want to remove the Wuhan Infected since they are unable to travel. 
    ChinaConfirmationData <- read.csv( "../../Data/ChinaCity_Confirmed_Map_0115_0318.csv" )
    wuhanInfected <- ChinaConfirmationData[match( "Wuhan", ChinaConfirmationData$City_EN ),
                                           match( "T_C_0201", names(ChinaConfirmationData)) - 1
                                           + as.numeric( substr(t0DateChina, 4,5 ) )]
        
    chinaPeriodOne <- chinaData[PeriodT0ColChina+Lag] - wuhanInfected
    chinaPeriodTwo <- chinaData[PeriodT1Col]


    
    italyPeriodOne <- italyData[PeriodT0ColItaly+Lag]
    italyPeriodTwo <- italyData[PeriodT1Col]

    spainPeriodOne <- spainData[PeriodT0ColItaly+Lag]
    spainPeriodTwo <- spainData[PeriodT1Col]
    germanPeriodOne <- germanData[PeriodT0ColItaly+Lag]
    germanPeriodTwo <- germanData[PeriodT1Col]
    ukPeriodOne <- ukData[PeriodT0ColItaly+Lag]
    ukPeriodTwo <- ukData[PeriodT1Col]

    chinaPop <- chinaData$Population..2020.
    italyPop <- italyData$Population..2020.
    spainPop <- spainData$Population..2020.
    germanyPop <- germanData$Population..2020.
    ukPop <- ukData$Population..2020.

    stateInfection$Confirmed <- stateInfection[,T1ColUS]
    
    MergedData <- merge( x = mDF[,c("Combined_Key","ChinaTravel","ItalyTravel","SpainTravel", "GermanyTravel", "UKTravel")], y = stateInfection[,c("Combined_Key","Confirmed","Pop")], by="Combined_Key" ) %>% filter( Confirmed > 0 )


    ## RHS: $ R_{c,0} M_i \int_{T_0}^{T_1} \frac{1}{N_c - R_{c,t} }$ without lags.
    ## RHS: $ R_{c,k} M_i \int_{T_0}^{T_1-k} \frac{1}{N_c - R_{c,t} }$ with lags. 

    ## Note that our M_i contains travel over the interval T_0 to T_1,
    ## so we sum from T_0 to T_1 - k, and divide by the total length
    ## from T_0 to T_1.
    reg <- function( sourceP1, Travel, sourceTSeriesCum, sourcePop ){
        exp( log( as.numeric(Travel) ) + log(as.numeric(sourceP1)))*sum( exp( - log( sourcePop - as.numeric(sourceTSeriesCum) ) ) )
    }
    

    ## chinaData and Italy data contains the cumulative infections:
    ## R_c, so sourceTSeriesCum is just the array of the cumulative
    ## from T_0 to T_1-k. But length is T_0 to T_1
    if( useItaly ){
        MergedData$Reg <- reg(  chinaPeriodOne, MergedData$ChinaTravel,
                              chinaData[PeriodT0ColChina:(PeriodT1Col-Lag)] - wuhanInfected,
                              chinaPop ) +
            reg( italyPeriodOne, MergedData$ItalyTravel,
                italyData[PeriodT0ColItaly:(PeriodT1Col-Lag)], italyPop)
        if( fullEU ){
            MergedData$Reg <- MergedData$Reg + 
                reg( spainPeriodOne, MergedData$SpainTravel,
                    spainData[PeriodT0ColItaly:(PeriodT1Col-Lag)], spainPop) +
                reg( germanPeriodOne, MergedData$GermanyTravel,
                    germanData[PeriodT0ColItaly:(PeriodT1Col-Lag)], germanyPop) +
                reg( ukPeriodOne, MergedData$UKTravel,
                    ukData[PeriodT0ColItaly:(PeriodT1Col-Lag)], ukPop)
        }
    }else{
        MergedData$Reg <- reg(  chinaPeriodOne, MergedData$ChinaTravel,
                              chinaData[PeriodT0ColChina:(PeriodT1Col-Lag)] - wuhanInfected, chinaPop )
    }
    
    ## For computing Reg for Iceland.
    iceEntryChina <- iceLandEntry$useMe[match("China", iceLandEntry$X2020)]
    iceEntryItaly <- iceLandEntry$useMe[match("Italy", iceLandEntry$X2020)]
    iceEntrySpain <- iceLandEntry$useMe[match("Spain", iceLandEntry$X2020)]
    iceEntryGermany <- iceLandEntry$useMe[match("Germany", iceLandEntry$X2020)]
    iceEntryUK <- iceLandEntry$useMe[match("United Kingdom", iceLandEntry$X2020)]

    # 67
    iceLandData <- worldInfection[match("Iceland", worldInfection$Country.Region),]

    icePop <- iceLandData$Population..2020.

    ## We use only the First wave of deCODE From Mar 15-19
    
    deCodePercent <- (48/5490)
    icePeriodTwo <- 0.0
    if( useDECODE ){
        icePeriodTwo <- (deCodePercent*icePop / iceLandData$X3.15.20)*iceLandData[PeriodT1Col-Lag]
    }else{
        icePeriodTwo <- iceLandData[PeriodT1Col-Lag]*2
    }

    if( useItaly ){
        denom <- reg(  chinaPeriodOne, iceEntryChina,
                              chinaData[PeriodT0ColChina:(PeriodT1Col-Lag)] - wuhanInfected,
                              chinaPop ) +
            reg( italyPeriodOne, iceEntryItaly,
                italyData[PeriodT0ColItaly:(PeriodT1Col-Lag)], italyPop)
        if( fullEU ){
            denom <- denom + 
                reg( spainPeriodOne, iceEntrySpain,
                    spainData[PeriodT0ColItaly:(PeriodT1Col-Lag)], spainPop) +
                reg( germanPeriodOne, iceEntryGermany,
                    germanData[PeriodT0ColItaly:(PeriodT1Col-Lag)], germanyPop) +
                reg( ukPeriodOne, iceEntryUK,
                    ukData[PeriodT0ColItaly:(PeriodT1Col-Lag)], ukPop)
        }
        
        iceAlphaBetaGamma <- icePeriodTwo / (denom)
    }else{

        iceAlphaBetaGamma <- icePeriodTwo / (reg( chinaPeriodOne, iceEntryChina, chinaData[PeriodT0ColChina:(PeriodT1Col-Lag)] - wuhanInfected, chinaPop ))
    }

    MergedData$alpha <- MergedData$Confirmed / (MergedData$Reg * as.numeric(iceAlphaBetaGamma) )

    MergedData$AlphaC <- 1.0 - MergedData$alpha##format( 1.0 - as.numeric(MergedData$alpha), digits=3)

    MergedData$InfectionsPerConfirmed <- MergedData$AlphaC / MergedData$alpha
        ##format( as.numeric(MergedData$AlphaC) / as.numeric(MergedData$alpha), digits=3)

    MergedData$iceAlphaBetaGamma <- as.numeric(iceAlphaBetaGamma)
    
    MergedData$Reg <- MergedData$Reg * as.numeric(iceAlphaBetaGamma)

    MergedData
}

mDF <- read.csv( "../../Data/countyLevelTravelData_US.csv", sep=";")
mDF2 <- read.csv("../../Data/DataFinal/EUtravel.csv" )[1:30,]
mDF2$SpainSum <- mDF2$SpainJan + mDF2$SpainFeb
mDF2$GermanySum <- mDF2$GermanyJan + mDF2$GermanyFeb
mDF2$UKSum <- mDF2$UKJan + mDF2$UKFeb
mDF2$Combined_Key <- paste( mDF$County, "US", sep=", ")

statePop <- read.csv("../../Data/US_CountyPop_Top30Ports.csv", sep=";")
mDF$Combined_Key <- paste( mDF$County, "US", sep=", ")
statePop$Combined_Key <- paste( statePop$County, "US", sep=", ")

mDF <- merge( mDF[,c("Combined_Key","ChinaSum","ItalySum")], mDF2[,c("Combined_Key", "SpainSum", "GermanySum", "UKSum" )], by="Combined_Key" )
## 60 is the number of days in the travel data.
mDF$ChinaTravel <- mDF$ChinaSum/60.0
mDF$ItalyTravel <- mDF$ItalySum/60.0
mDF$SpainTravel <- mDF$SpainSum/60.0
mDF$GermanyTravel <- mDF$GermanySum/60.0
mDF$UKTravel <- mDF$UKSum/60.0
iceLandEntry$useMe <- (iceLandEntry$Jan+iceLandEntry$Feb)/60.0

t0DateChina <- "02-23-2020"
t0DateItaly <- "02-23-2020"



## Table 2 - Lag 5
cutoffDate <- "03-10-2020"
Lag <- 5
## For Reference: Here are the Three Booleans in EstimateAlpha()
## excludeEarly, useItaly, useDECODE, Lag, fullEU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Only China
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, FALSE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)

## Full EU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Table 2: Lag 8
cutoffDate <- "03-13-2020"
Lag <- 8


## excludeEarly, useItaly, useDECODE, Lag, fullEU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Only China
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, FALSE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)

## Full EU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)

## Reproduction of the Figure
cutoffDate <- "03-10-2020"
Lag <- 5


## excludeEarly, useItaly, useDECODE )
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)

## This will need to be edited as the data becomes more riched, and
## hopefully automated.  For now, this just makes the county names on
## the plot not seem to crowded.- If you want to look at King County,
## it needs to be added to this file, or change label=shortnames to
## label=Combined_Key

## For Lag = 5
MergedData$shortnames <- c("Broward, FL", "Clark, NV", "Cook, IL", "Fulton, GA", "Harris, TX", "Hillsborough, FL", "Honolulu, HI", "Los Angeles, CA", "Maricopa AZ", "New York City, NY", "Philadelphia, PA", "Ramsey, MN", "San Diego, CA", "San Francisco, CA", "Suffolk, MA" )


## Our standard error estimates
## sqrt(vcov(alphaReg)[1,1] )

## If you wish to look at each city individually
## stargazer( MergedData[,c(1,4,7,8,9)], summary=FALSE )

temp <- function(x) {alphaReg$coef[1]*x}

pdf( "ScatterAlpha_Full.pdf")
fit <- lm(Confirmed ~ Reg - 1, data = MergedData)
ggplot(MergedData, aes(x=Reg, y=Confirmed)) +
    geom_point(shape=1) +

      # Use hollow circles
    stat_function(fun = temp,
                  color = "red", size=.5) +

     geom_label_repel(aes(label = shortnames),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') + 
    ylab("Confirmed Cases") +
    xlab("Estimated Infectious")
dev.off()



## Reproduction of the second Figure
cutoffDate <- "03-13-2020"
Lag <- 8

t0DateChina <- "02-23-2020"
t0DateItaly <- "02-23-2020"
## excludeEarly, useItaly, useDECODE )
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)

## This will need to be edited as the data becomes more riched, and
## hopefully automated.  For now, this just makes the county names on
## the plot not seem to crowded.- If you want to look at King County,
## it needs to be added to this file, or change label=shortnames to
## label=Combined_Key



## For Lag = 8
MergedData$shortnames <- c("Broward, FL", "Clark, NV", "Cook, IL", "Dallas, TX", "Essex, NJ", "Fulton, GA", "Harris, TX", "Hillsborough, FL", "Honolulu, HI", "Los Angeles, CA", "Maricopa AZ", "Miami-Dade, FL", "Multnomah, OR", "New York City, NY", "Philadelphia, PA", "Ramsey, MN", "San Diego, CA", "San Francisco, CA", "Suffolk, MA", "Wayne, MI", "Whatcom, WA" )

## Our standard error estimates
## sqrt(vcov(alphaReg)[1,1] )

## If you wish to look at each city individually
## stargazer( MergedData[,c(1,4,7,8,9)], summary=FALSE )

temp <- function(x) {alphaReg$coef[1]*x}

pdf( "ScatterAlpha_Full_Lag8.pdf")
fit <- lm(Confirmed ~ Reg - 1, data = MergedData)
ggplot(MergedData, aes(x=Reg, y=Confirmed)) +
    geom_point(shape=1) +

      # Use hollow circles
    stat_function(fun = temp,
                  color = "red", size=.5) +

     geom_label_repel(aes(label = shortnames),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') + 
    ylab("Confirmed Cases") +
    xlab("Estimated Infectious")
dev.off()





matFive <- matrix( 0, 19, 29)
Lag <- 5


## Note that because Iceland features its first infection on Feb 28,
## Lagging with 5 days means we can only really start at Mar 6 before
## we start getting some nonsense results.
for( i in 6:19 ){
    ## Construct the string giving the date that will be used to find the file.
    if(nchar(  as.character(i) ) == 1){
        cutoffDate <- paste( "03-0", as.character(i),"-2020", sep = "")
    }else{
        cutoffDate <- paste( "03-", as.character(i),"-2020", sep = "")
    }
    
    

    for( j in 1:29){
        if(nchar(  as.character(j) ) == 1){
            t0DateChina <- paste( "02-0", as.character(j),"-2020", sep = "")
        }else{
            t0DateChina <- paste( "02-", as.character(j),"-2020", sep = "")
        }

        MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE )

        alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)

        matFive[i,j] <- format(  alphaReg$coef[1], digits=4) ##format( alphaReg$coef[1], digits=4)
    }
}

matEight <- matrix( 0, 19, 29)
Lag <- 8


## Note that because Iceland features its first infection on Feb 28,
## Lagging with 5 days means we can only really start at Mar 6 before
## we start getting some nonsense results.
for( i in 9:19 ){
    ## Construct the string giving the date that will be used to find the file.
    if(nchar(  as.character(i) ) == 1){
        cutoffDate <- paste( "03-0", as.character(i),"-2020", sep = "")
    }else{
        cutoffDate <- paste( "03-", as.character(i),"-2020", sep = "")
    }
    
    

    for( j in 1:29){
        if(nchar(  as.character(j) ) == 1){
            t0DateChina <- paste( "02-0", as.character(j),"-2020", sep = "")
        }else{
            t0DateChina <- paste( "02-", as.character(j),"-2020", sep = "")
        }
        

        MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE )

        alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)

        matEight[i,j] <- format(  alphaReg$coef[1], digits=4) ##format( alphaReg$coef[1], digits=4)
    }
}


## I apply a regular expression replace on this table to produce the
## date, but this is contentwise the same.
stargazer( cbind( 6:19, matFive[6:19,c(1,5,10,15,20,25,29)] ) )
stargazer( cbind( 9:19, matEight[9:19,c(1,5,10,15,20,25,29)] ) )

## Vec stores lags based around Mar 5 - Same as Mar 10 with 5 day lag. 
vec <- numeric(16)
start <- 5
for( l in 0:15 ){
    ## Construct the string giving the date that will be used to find the file.
    if(nchar(  as.character(start+l) ) == 1){
        cutoffDate <- paste( "03-0", as.character(start+l),"-2020", sep = "")
    }else{
        cutoffDate <- paste( "03-", as.character(start+l),"-2020", sep = "")
    }
    Lag <- l

    t0DateChina <- "02-23-20"
    

    MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE )

    alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)

    vec[l+1] <- format(  alphaReg$coef[1], digits=3) ##format( alphaReg$coef[1], digits=4)
}

stargazer( cbind( 0:15, vec ) )





## Now we get the results including Seattle (King County)
## Lag 5
cutoffDate <- "03-10-2020"
Lag <- 5
## For Reference: Here are the Three Booleans in EstimateAlpha()
## excludeEarly, useItaly, useDECODE, Lag, fullEU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, FALSE, TRUE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Only China
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, FALSE, FALSE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)

## Full EU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, FALSE, TRUE, TRUE, Lag, TRUE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Table 2: Lag 8
cutoffDate <- "03-13-2020"
Lag <- 8


## excludeEarly, useItaly, useDECODE, Lag, fullEU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, FALSE, TRUE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Only China
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, FALSE, FALSE, TRUE, Lag, FALSE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)

## Full EU
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, FALSE, TRUE, TRUE, Lag, TRUE  )

alphaReg <- lm( Confirmed ~ Reg - 1, data =MergedData)
alpha <- alphaReg$coef[1]

print(alpha)
print( (1-alpha)/alpha)


## Reproduction of "Table x: Estimated Total Infected By County"

stateInfection <- read.csv("../../Data/time_series_covid19_confirmed_US.csv" )
stateInfection <- merge( x=statePop[,c("Pop", "Combined_Key")], y = stateInfection,
                        by="Combined_Key" )

cutoffDate <- "03-10-2020"
Lag <- 5


t0DateChina <- "02-23-2020"
t0DateItaly <- "02-23-2020"
## excludeEarly, useItaly, useDECODE )
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE  )

MergedData$shortnames <- c("Broward, FL", "Clark, NV", "Cook, IL", "Fulton, GA", "Harris, TX", "Hillsborough, FL", "Honolulu, HI", "Los Angeles, CA", "Maricopa AZ", "New York City, NY", "Philadelphia, PA", "Ramsey, MN", "San Diego, CA", "San Francisco, CA", "Suffolk, MA" )
Mar15Date <- match( "X3.15.2020", names(stateInfection ) )
Mar20Date <- match( "X3.20.2020", names(stateInfection ) )
stateInfection <- read.csv("../../Data/time_series_covid19_confirmed_US.csv" ) %>%
    filter( !(Combined_Key %in% c( "King, Washington, US" ) ) )

MergedData <- merge( MergedData, stateInfection, by="Combined_Key" )
           
Mar15Date <- match( "X3.15.2020", names(MergedData ) )
Mar20Date <- match( "X3.20.2020", names(MergedData ) )

## Haris Texas has a 0 element on the cumulative for Mar 20.
## So we impute the Mar 19 number to be safe.
MergedData[match("Harris, Texas, US",MergedData$Combined_Key),Mar20Date] <-
    MergedData[match("Harris, Texas, US",MergedData$Combined_Key),Mar20Date-1]

dat <- sort( matFive[6:19,] )
minNum <- as.numeric(dat[length(dat)*.9])
meanNum <- mean( as.numeric(dat))
maxNum <- as.numeric(dat[length(dat)*.1])

MergedData$Mar15Min <- round( MergedData[,Mar15Date] / minNum)
MergedData$Mar15Mean <- round( MergedData[,Mar15Date] / meanNum)
MergedData$Mar15Max <- round( MergedData[,Mar15Date] / maxNum)

MergedData$Mar20Min <- round( MergedData[,Mar20Date] / minNum)
MergedData$Mar20Mean <- round( MergedData[,Mar20Date] / meanNum)
MergedData$Mar20Max <- round( MergedData[,Mar20Date] / maxNum)

stargazer( MergedData[,c("shortnames", "X3.15.2020", "Mar15Min", "Mar15Mean", "Mar15Max", "X3.20.2020", "Mar20Min", "Mar20Mean", "Mar20Max")], summary = FALSE )


cutoffDate <- "03-13-2020"
Lag <- 8


t0DateChina <- "02-23-2020"
t0DateItaly <- "02-23-2020"
## excludeEarly, useItaly, useDECODE )
MergedData <- EstimateAlpha( mDF, statePop, worldInfection, t0DateChina, t0DateItaly, cutoffDate, iceLandEntry, TRUE, TRUE, TRUE, Lag, TRUE  )

MergedData$shortnames <- c("Broward, FL", "Clark, NV", "Cook, IL", "Dallas, TX", "Essex, NJ", "Fulton, GA", "Harris, TX", "Hillsborough, FL", "Honolulu, HI", "Los Angeles, CA", "Maricopa AZ", "Miami-Dade, FL", "Multnomah, OR", "New York City, NY", "Philadelphia, PA", "Ramsey, MN", "San Diego, CA", "San Francisco, CA", "Suffolk, MA", "Wayne, MI", "Whatcom, WA" )
Mar15Date <- match( "X3.15.2020", names(stateInfection ) )
Mar20Date <- match( "X3.20.2020", names(stateInfection ) )
stateInfection <- read.csv("../../Data/time_series_covid19_confirmed_US.csv" ) %>%
    filter( !(Combined_Key %in% c( "King, Washington, US" ) ) )

MergedData <- merge( MergedData, stateInfection, by="Combined_Key" )
           
Mar15Date <- match( "X3.15.2020", names(MergedData ) )
Mar20Date <- match( "X3.20.2020", names(MergedData ) )

## Haris Texas has a 0 element on the cumulative for Mar 20.
## So we impute the Mar 19 number to be safe.
MergedData[match("Harris, Texas, US",MergedData$Combined_Key),Mar20Date] <-
    MergedData[match("Harris, Texas, US",MergedData$Combined_Key),Mar20Date-1]

MergedData[match("Miami-Dade, Florida, US",MergedData$Combined_Key),Mar20Date] <-
    mean(c(MergedData$X3.18.2020[match("Miami-Dade, Florida, US",MergedData$Combined_Key)],
           MergedData$X3.21.2020[match("Miami-Dade, Florida, US",MergedData$Combined_Key)]))

dat <- sort( matEight[9:19,] )
minNum <- as.numeric(dat[length(dat)*.9])
meanNum <- mean( as.numeric(dat))
maxNum <- as.numeric(dat[length(dat)*.1])

MergedData$Mar15Min <- round( MergedData[,Mar15Date] / minNum)
MergedData$Mar15Mean <- round( MergedData[,Mar15Date] / meanNum)
MergedData$Mar15Max <- round( MergedData[,Mar15Date] / maxNum)

MergedData$Mar20Min <- round( MergedData[,Mar20Date] / minNum)
MergedData$Mar20Mean <- round( MergedData[,Mar20Date] / meanNum)
MergedData$Mar20Max <- round( MergedData[,Mar20Date] / maxNum)

stargazer( MergedData[,c("shortnames", "X3.15.2020", "Mar15Min", "Mar15Mean", "Mar15Max", "X3.20.2020", "Mar20Min", "Mar20Mean", "Mar20Max")], summary = FALSE )

MergedData$PercentageMar15 <- MergedData$Mar15Max / MergedData$Pop
MergedData$PercentageMar20 <- MergedData$Mar20Max / MergedData$Pop

