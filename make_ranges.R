
#this function finds contiguous elements in a vector
splitter = function(diffs){
  pos = list(-1)
  neg = list(-1)
  Vec = diffs
  Breaks <- c(0, which(diff(Vec) != 1), length(Vec))
  if(length(Breaks)==2){
    res = list(Vec[(Breaks[1] + 1):Breaks[2]])
  }else{
    res = sapply(seq(length(Breaks) - 1), function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
  }
  
  return(res)
}


#this function uses the above function to find contiguous elements in a given list of dm elements. This will be called on
#LL@Intervals$pos, LL@Intervals$neg. l is the value of lambda corresponding to the dmrs. It comes from:
#LL@lambda[[i]]. chr is the chrosome of the data. It returns an IRanges object
make.ranges.single = function(pos, neg,l, chr){
  require(IRanges)
  
  DMRS=IRanges()
  for(i in 1:length(pos)){
    startsp = pos[[i]][1]
    stopsp = pos[[i]][length(pos[[i]])]
    
    Irang2pos = IRanges(startsp, stopsp)
    elementMetadata(Irang2pos)= list("lambda"=rep(l, length(1)), "direction"=rep("pos",length(1) ), "chr"=rep(chr,length(1) ))
    
    DMRS=unlist(IRangesList(DMRS, Irang2pos))
  }
  
  if(length(neg)){
    for(i in 1:length(neg)){
      
      startsp = neg[[i]][1]
      stopsp = neg[[i]][length(neg[[i]])]
      
      Irang2neg = IRanges(startsp, stopsp)
      elementMetadata(Irang2neg)= list("lambda"=rep(l, length(1)), "direction"=rep("neg",length(1) ), "chr"=rep(chr,length(1) ))
      
      DMRS=unlist(IRangesList(DMRS, Irang2neg))
    }
  }
  
  return(DMRS)
  
}

###this is a wrapper function that loops through the intervals and creates the GRanges objects

make.ranges = function(multi, chr){
  pos = multi@Intervals$pos
  neg = multi@Intervals$negs
  ls = multi@lambdas
  DMRS= IRanges()
  for(i in 1:length(pos)){
    ranges=make.ranges.single(pos[[i]], neg[[i]], ls[i],chr)
    DMRS=unlist(IRangesList(DMRS,ranges))
  }
  
  DMRS = GRanges(mcols(DMRS)$chr, IRanges(DMRS@start, DMRS@start+DMRS@width), strand="+", mcols(DMRS)[,c(1,2)])
  return(DMRS)
}

