#Mean allelic richness estimated by rarefaction

AllelicRichRaref <- function(data){
  V = c()
  vector = c()
  maxsample <- nrow(data)
  for(i in 1:maxsample){                                                  #Resampling
    G <- 2*i                                                              #Number of genes
    loci <- ncol(data)/2
    for(j in 1:loci){                                                     #Analysis each locus
      data2 <- data[ ,((j*2)-1):(j*2)]                                    #Select the columns of the locus
      K <- as.data.frame(table(unlist(data2)))                            #Calculate the frequency of each allele
      K$Var1 <- as.numeric(K$Var1)
      my_vector = c()
      for(k in K$Var1){                                                   #Analysis each allele                                
        Ni <- nrow(data)*2
        Nik <- K[K$Var1==k,]$Freq
        numerator <- choose((Ni-Nik),G)                                   #Number of combinations of G genes that not include the allele k
        denominator <- choose((Ni),G)                                     #Number of possible combinations of G genes taken from the sample of Ni genes
        Pi <- numerator/denominator                                       #Probability that allele k does not occur in a sample of G genes chosen as reference
        nPi <- 1-Pi                                                       #Probability that allele k occurs
        my_vector[length(my_vector)+1] <- nPi
      }
      E <- sum(my_vector)                                                 #Allelic richness at one locus
      vector[length(vector)+1] <- E
    }
    A <- mean(vector)                                                     #Mean of the allelic richness
    V[length(V)+1] <- A                                                   #List of the means of each resample
  }
  allelicrich <- V
  table <- data.frame(allelicrich)
  table$SampleSize <-c(1:maxsample)
  table[,c(2,1)]
}
