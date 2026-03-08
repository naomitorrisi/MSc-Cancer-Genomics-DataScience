# Functions required for Lab_2: Statistical analysis of quantitative proteomics data
# Author: Pedro R. Cutillas



ksear.s <- function(df.fold, ks_db){

  # df.fold == dataset of fold changes or
  # ks_db == database of kinase-substrate relationships
  #       possibilities are "signor" "pdts", "pSite"
  df.fold[,1] <- gsub("..",");",df.fold[,1],fixed = T)
  df.fold[,1] <- gsub(".","(",df.fold[,1],fixed = T)
  values.all <- na.omit(as.numeric(subset(df.fold[,2], df.fold[,2]!=0)))
  nc <- ncol(df.fold)
  df.ks <- get.ks.set(ks_db)
  nr <- nrow(df.ks)
  zscores <- numeric(nr)
  pvalues <- numeric(nr)
  msites <- numeric(nr)
  kinases <- character(nr)
  allsites <- character(nr)
  r=1
  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,3])
        ss <- c(unlist(strsplit(substrates,";")))
        start.time <- Sys.time()
        df.xx <- subset(df.fold,df.fold[,1] %in% paste(ss,";",sep=""))
        sites <-paste(unlist(df.xx[,1]),collapse=";")
        if (nrow(df.xx)>2){
         # df.xx$prot.group <- kinase
          #sites.x <- data.frame(site=df.xx[,1])
          #sites.x$kinase <- kinase
          #all.sites <- rbind(all.sites,sites.x)

            myvalues <- na.omit(as.numeric(subset(df.xx[,2],df.xx[,2]!=0)))
            pval <- 1
            tryCatch({
              myks <- ks.test(values.all,myvalues,exact=F)
              pval <- myks$p.value
            }, error=function(e){}
            )
            mysd <- sd(values.all,na.rm = T)
            mymedian <- median(myvalues, na.rm = T)
            mymedian.all <- median(values.all,na.rm = T)
            sd.all <- sd(values.all)
            msites[r] <- length(myvalues)
            zscores[r] <- ((mymedian-mymedian.all)*((sqrt(msites[r])))/mysd)
            pvalues[r] <- pval
            kinases[r] <- kinase
            allsites[r] <- sites
          }
        }
      }
  }

  df.out <- data.frame(kinases,zscores,pvalues,m=msites,allsites)
  df.out <- subset(df.out,df.out$m>1)
  df.out$qval <- p.adjust(df.out$pvalues,method = "BH")
  return(df.out)

  ##################################
}


compare.by.limma <- function(df.to.compare, control.samples, test.samples){


  # Compare by limma
  #
  # first column is protein or ppsite names
  # first set of samples are control
  # second set of samples are the test
  library(limma)
  # df.to.compare <- df.ppindex[,1:9]
  nc <- ncol(df.to.compare)
  #replicates <- (nc-1)/2
  control.samples <- intersect(control.samples,colnames(df.to.compare))
  test.samples <- intersect(test.samples,colnames(df.to.compare))
  df.s <- df.to.compare[,c(control.samples,test.samples)]
  df.s1 <- data.frame(outcome=matrix(nrow=length(control.samples)))
  df.s2 <- data.frame(outcome=matrix(nrow=length(test.samples)))
  df.s1$outcome <- "control"
  df.s2$outcome <- "test"
  df.ss <- rbind(df.s1,df.s2)
  des <- factor(ifelse(df.ss$outcome=="control" ,"1",
                       "2"))
  facna <- addNA(des)
  design <- model.matrix(~ 0+factor(c(facna)))
  colnames(design) <- c("control","test")
  contrast.matrix <- makeContrasts(test-control,
                                   levels=design)
  #df.s ==== proteomics data
  fit <- lmFit(df.s,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  pvals <- data.frame(fit2$p.value)
  fvals <- data.frame(fit2$coefficients)
  df.xx <- data.frame(protein=df.to.compare[,1],
                      differnces=fvals,
                      pvalues=pvals)
  colnames(df.xx) <- c("protein","difference.test.vs.control","pvalues")

  df.xx$FDR <- p.adjust(df.xx$pvalues, method = "fdr")

  return(df.xx)

}


kinase.substrate.enrichment <- function(dfx, ks_db,is.ksea=TRUE){



  # df.fold == dataset of fold changes or
  # ks_db == database of kinase-substrate relationships
  #       possibilities are "edges", "ctams" "pdts", "pSite"
  #
  # returns 6 data frames:
  #         results.zscores,
  #         results.distance
  #         results.pvalues,
  #         results.m,
  #         results.q
  #         results.sites
  dfx <- data.frame(dfx)
  rownames(dfx) <- make.names( dfx[,1], unique = T)
  if (is.ksea==TRUE){
    column.with.priors <- 3
  }else{
    #Convert phosphopeptides to genes for phosphoproteomics data when matching to proteomics sets
    x1 <-  grepl("(", as.character(dfx[,1]),fixed=T)
    if (TRUE %in% x1){ # True if phosphoproteomics data, false if proteomics data
      dfx[,1] <- unlist(lapply(dfx[,1],
                               function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1]))
      column.with.priors <- "genes"
    }
  }

  dfx[,1] <- gsub("..",");",dfx[,1],fixed = T)
  dfx[,1] <- gsub(".","(",dfx[,1],fixed = T)
  nc <- ncol(dfx)


  ###################################################
  # Alternatives to get the ks_db:
  df.ks <- protools2::protein_and_ks_sets[[ks_db]]

  #df.ks <- read.csv("pdts") # protools2::protein_and_ks_sets[[ks_db]]

  ##############################

  nr <- nrow(df.ks)
  results.pvalues <-numeric(nr)
  results.zscores <-numeric(nr)
  results.q <-numeric(nr)
  results.m <-numeric(nr)
  results.distance <-numeric(nr)
  results.sites <-character(nr)
  kinases <- character(nr)
  r=1
  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,column.with.priors])
        ss <- c(unlist(strsplit(substrates,";")))
        start.time <- Sys.time()

        df.xx <- subset(dfx,dfx[,1] %in%  paste(ss,";",sep=""))

        sites <-paste(unlist(rownames(df.xx)),collapse=";")


        if (nrow(df.xx)>2){
          df.xx$prot.group <- kinase
          sites.x <- data.frame(site=df.xx[,1])
          sites.x$kinase <- kinase
          #all.sites <- rbind(all.sites,sites.x)
          c=2
          ds <- numeric(nc-1)
          pvals <- numeric(nc-1)
          zscores <- numeric(nc-1)
          ms<- numeric(nc-1)
          qs <- numeric(nc-1)

          values.all <- na.omit(as.numeric(subset(dfx[,c], dfx[,c]!=0)))
          myvalues <- na.omit(as.numeric(subset(df.xx[,c],df.xx[,c]!=0)))
          pval <- 1
          tryCatch({
            myks <- ks.test(values.all,myvalues)
            pval <- myks$p.value
          }, error=function(e){}
          )

          m <- 0
          q <- 0
          m <- nrow(df.xx)
          mysd <- sd(values.all,na.rm = T)
          mymedian <- median(myvalues, na.rm = T)
          mymedian.all <- median(values.all,na.rm = T)
          sd.all <- sd(values.all)

          results.distance[r] <- mymedian-mymedian.all
          results.pvalues[r] <- pval
          results.zscores[r] <- ((mymedian-mymedian.all)*((sqrt(m)))/mysd)
          results.m[r] <- m
          # results.q[r] <- qs
          results.sites[r] <- as.character(sites)
          kinases[r] <- kinase

        }
      }
    }
  }

  xx <- data.frame(kinases,
                   zscores=results.zscores,
                   pvalues=results.pvalues,
                   m=results.m,
                   #q=results.q,
                   distance=results.distance,
                   sites=results.sites,
                   kinase_dataset =paste0(kinases,"_",ks_db))

  xx <- subset(xx,xx$m>1)

  #xx <- xx[order(xx$zcores),]
  xx$qvalue <- p.adjust(xx$pvalues,method = "fdr")
  #head(xx[order(xx$pvalues),],n=20)
  return(xx)
  ##################################
}


