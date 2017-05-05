###function that takes data from bw objects containing both samples (normal type and cancer type) and creates a data frame
###for the interval finder.

###file.name.stem is the name of the file with the "stem cell" data (this would be normal iff comparing normal v cancer)
###file.name.can is the name of the file with the cancer data in i
###chr gives the chromosome to be investigated. The rest of the data is thrown away. chr can also be set to "all" which
###keeps the entire genome

prep.data = function(file.name.stem,file.name.can, chr, file.type="bw", convert=TRUE){
  require(rtracklayer)
  
  ###import the data. File type defines what kind of file the input comes in. .txt corresponds to methyl data from TCGA
  if(file.type=="bw"){
    D_stem = import(file.name.stem, format="bw")
    if(chr!="all"){
      D_stem_chr = D_stem[which(as.vector(D_stem@seqnames)==chr),]
    }else{
      D_stem_chr = D_stem
    }
    
    ###code assumes score contains the the methylation and start contains the loci corresponding to the methylation value
    #has the methylation score
    stem_turn_chr=D_stem_chr@elementMetadata@listData$score
    #has the loci
    stem_x = D_stem_chr@ranges@start
    
    
    D_can = import(file.name.can, format="bw")
    if(chr!="all"){
      D_can_chr = D_can[which(as.vector(D_can@seqnames)==chr),]
    }else{
      D_can_chr = D_can
    }
    can_turn_chr=D_can_chr@elementMetadata@listData$score
    can_x = D_can_chr@ranges@start
  }else if(file.type=='txt'){
    D_stem=read.delim(file.name.stem,sep="\t",stringsAsFactors=FALSE)
    rms=which(D_stem$Start<0)
    
    D_can=read.delim(file.name.can,sep="\t",stringsAsFactors=FALSE)
    rms1 = which(D_can$Start<0)
    
    rms2n = which(is.na(D_stem$Beta_value))
    rms2c = which(is.na(D_can$Beta_value))
    
    D_can = D_can[-c(rms, rms2n, rms2c),]
    D_can=D_can[order(D_can$Start),]
    
    D_stem=D_stem[-c(rms, rms2n, rms2c),]
    D_stem= D_stem[order(D_stem$Start),]
    
    D_can_chr=D_can[which(D_can$Chromosome==chr),]
    D_stem_chr=D_stem[which(D_stem$Chromosome==chr),]
    
    can_x = D_can$Start
    stem_x = D_stem$Start
    
    can_turn_chr = D_can_chr$Beta_value
    stem_turn_chr= D_stem_chr$Beta_value
    
    if(convert=="TRUE"){
      Beta_c = can_turn_chr
      Beta_n= stem_turn_chr
      M_c=log2(Beta_c/(1-Beta_c))
      M_n= log2(Beta_n/(1-Beta_n))
      can_turn_chr=M_c
      stem_turn_chr=M_n
    }
  }
  
  ###this finds all of the loci which are contained in both data sets. Not actually necessary to do this but it makes
  ###things slightly easier at the estimation stage
  shared_x = which(can_x%in%stem_x)
  can_x_full = can_x
  can_x = can_x_full[shared_x]
  
  shared_x_stem = which(stem_x%in%can_x)
  stem_x_full = stem_x
  stem_x = stem_x[shared_x_stem]
  
  ###stem x now contains the shared loci, so we subset the methylation values to correspond with the shared loci
  can_turn_chr_full = can_turn_chr
  can_turn_chr=can_turn_chr[shared_x]
  stem_turn_chr_full = stem_turn_chr
  stem_turn_chr = stem_turn_chr[shared_x_stem]
  
  ###assigning names that correspond with names from the interval estimation function
  x.n = stem_x
  x.c = x.n
  y.n= stem_turn_chr
  y.c = can_turn_chr
  
  return(list("x.n"=x.n, "x.c"=x.c, "y.n"=y.n, "y.c"=y.c))
}
