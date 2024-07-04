# The code is written by Kwanghun Choi (2023-09-26)
## This is the collection of codes from Kwon and Chung (2004)
### This collection is proto type and I will split functions into modules later


#### 1. Tree height estimation function
HGHT  <- 
    function(DBH)
    {
        result <- 32.818 * exp( - 24.455 / (DBH + 7.836) )
        return(result)
    }

#### 2. Crown ratio function
CR    <-   
    function(DBH, BA)
    {
        result <- 0.1733 * ( 1 / (1 + 0.0312 * BA) ) + 0.4414 * ( 1 - exp( -0.1864 * DBH ) )
        return(result)
    }

#### 3. Potential diameter growth function
PG    <- 
    function(DBH, SI, CR)
    {
        result <- 0.6114 - 0.0006 * DBH^(1.7820) + 0.0092 * SI * CR * DBH^(0.6384)  
        return(result)
    }

#### 4. Modifier function
MOD   <- 
    function(DBH, BA, BA_max = 90)
    {
        AD     <- mean(DBH, na.rm = TRUE)
        f_x    <- 4.1410 * ( 1 - exp(-0.1754 * DBH / AD) )^(2.9996) + 0.0042
        g_x    <- 3.5440 * ( AD + 1 )^(0.9998)

        result <- 1 - exp(- ( f_x * g_x ) * sqrt( max(( BA_max - BA ), 0) / BA ) )
        return(result)
    }

#### 5. Mortality function (I used MOR rather than M)
MOR   <- 
    function(DBH, PG, MOD)
    {
        result <- 1 / ( 1 + exp( 1.5238 + 5.2536 * ( PG * MOD )^(2.3811) + 0.0011 * DBH ) ) + 0.0210
        return(result)
    }

#### 6. Taper equation (h 는 벡테라는 점을 고려하여 다시 작성)
STR   <-
    function(DBH, H, interval=0.2, PARAM, SPP) # interval means height interval from the ground
    {
    # 1 = "Pinus koraensis" : Need tables for parameter sets (PARAM) of each species
    # PRM[1:3] = a0 - a2, PRM[4:8] = b1 - b5, PRM[9]  p, PRM[10] = FI
        # This equation is from Kozak 1988 model [Kozak, A. (2004) My last words on taper equations. The Forestry Chronicle., 80(4) 507-515]
        h           <- seq(0, H, by = interval) # vector with length(h)
        #PRM         <- PARAM |> filter(Species == SPP) |> subset(select(c(1:10))) |> t() |> c() # vector with length 10
        PRM         <- PARAM[,SPP] 
        Z           <- h / H # vector with length(h)
        # p           <- HI / H # number From parameter table 
        X           <- ( 1 - sqrt( h / H ) ) / ( 1 - sqrt( PRM[9] ) ) # vector with length(h)
        EXPONENT    <- PRM[4] * Z^2 + PRM[5] * log( Z + 0.001 ) + PRM[6] * sqrt(Z) + PRM[7] * exp(Z) +PRM[5] * ( DBH / H ) # vector with length(h)

        d_Kozak     <- PRM[1] * DBH^(PRM[2]) * PRM[3]^(DBH) * X^( EXPONENT ) # vector with length(h) # d is a vector from 1 to n for each h = 0 to ~1.
        dbh_Kozak   <- d_Kozak[c(round(1.2 / interval, 1) + 1)]

        # Calculate Volume following the Smalian's formula
        A_Kozak     <- pi/4 * d_Kozak^2
        A_Kozak_lag <- c(0, A_Kozak)[-c(length(A_Kozak)+1)]
        V_Somalian  <- sum(0.5*(A_Kozak+A_Kozak_lag) * interval)  + 0.5 * tail(A_Kozak, 1) * (H %% interval)
        return(data.frame(dbh_Kozak, V_Somalian))
    }
# 흉고직경:1.2m interval == 0.2일 때 7번째 d[7] 즉 round(1.2 / interval, 1) + 1 인덱스가 1부터 시작하므

# 임분별이기 때문에 임분 읽고 그 속에 유니크한 수종 뽑고 연도별 수종별 식 적
ForestPredictor <- 
    function(InputData, PARAM, TargetSite, InitialYear, SiteArea = 0.06, Period = 100, BA_max = 90, HghtInterval=0.2)
    {
        require(dplyr)
        temp_data      <- filter(InputData, grepl(TargetSite, No.))
        init_str       <- temp_data[, c(c(1:8), which(grepl(InitialYear, names(temp_data))))] |> 
                        mutate(across(contains(InitialYear), \(x) as.numeric(x)))
        temp_dbh       <- as.numeric(init_str[, grepl("직경", names(init_str))])
        temp_BA        <- sum(pi * ( temp_dbh / 2)^2, na.rm=T) / 10000 / SiteArea
        temp_height    <- as.numeric(init_str[, grepl("수고", names(init_str))])
        temp_spp       <- init_str[, "종명"]
        tree_str       <- do.call(bind_rows, apply(data.frame(temp_dbh, temp_height, temp_spp), 1, function(x){STR(DBH = as.numeric(x[1]), H = as.numeric(x[2]), SPP = x[3], PARAM = PARAM, interval = HghtInterval)}))
        tree_vol       <- tree_str[,2]
        result         <- bind_cols(temp_dbh, temp_height, tree_vol) |> setNames(paste0(c("직경", "수고", "재적"), formatC(as.numeric(InitialYear), width=3, flag="0")))

        for( i in c(1:Period))
        {
            
            temp_CR     <- CR(DBH = temp_dbh, BA=temp_BA)
            temp_PG     <- PG(DBH = temp_dbh, SI = 16, temp_CR)
            temp_MOD    <- MOD(DBH = temp_dbh, BA = temp_BA, BA_max = BA_max)
            temp_GR     <- temp_PG * temp_MOD
            temp_MOR    <- MOR(DBH = temp_dbh, PG = temp_PG, MOD = temp_MOD)

            temp_dbh    <- temp_dbh + temp_PG * temp_MOD
            temp_BA     <- sum(pi * ( temp_dbh / 2)^2, na.rm=T) / 10000 / SiteArea
            temp_height <- HGHT(DBH = temp_dbh)

            tree_str    <- do.call(bind_rows, apply(data.frame(temp_dbh, temp_height, temp_spp), 1, function(x){STR(DBH = as.numeric(x[1]), H = as.numeric(x[2]), SPP = x[3], PARAM = PARAM, interval = HghtInterval)}))
            tree_vol    <- tree_str[,2]
            temp_result <- bind_cols(temp_dbh, temp_height, tree_vol) |> setNames(paste0(c("직경", "수고", "재적"), formatC(as.numeric(InitialYear)+i, width=3, flag="0")))
            result      <- bind_cols(result, temp_result)
        }
        return(result)
    }


##### Data structure
# 임소반 정보, 수, 직경, 수고, 직경, 수고, 직경, 수고
# 에시 파일
InputExample <- read.csv("C:/R_storage/현장조사자료_재적건중량계산.csv")
ParamExample <- data.frame(잣나무 = c(0.9417, 1.0201, 0.997, 0.7667, -0.1556, 1.1785, -0.5214, 0.1235, 0.18, 0.980))
TargetSite   <- "A-Ⅰ"


Test <- ForestPredictor(InputData = InputExample, PARAM = ParamExample, TargetSite = TargetSite, InitialYear="06", Period = 100, BA_max = 90, HghtInterval = 0.2)
    

TestDT<- t((Test|>select(any_of(contains("직경"))))[1,])
plot(x=c(1:nrow(TestDT)), y=TestDT)

TestDT<- t((Test|>select(any_of(contains("수고"))))[1,])
plot(x=c(1:nrow(TestDT)), y=TestDT)

TestDT<- t((Test|>select(any_of(contains("재적"))))[1,])
plot(x=c(1:nrow(TestDT)), y=TestDT)

save.image("~/Downloads/ForestModel.rds")
