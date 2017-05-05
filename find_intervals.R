###these are the class definitions for the splines output which comes out of the the find_intervals function
multiscaless=setClass("multiscaless", slots=c(tau="numeric",Intervals="list",lambdas="numeric", splines="list", number="numeric"))

setMethod("show", "multiscaless", function(object) writeLines(paste(c("Multiscale Smoothing Spline object:","\n",
                                                                      "Spline Call:", paste(as.character(formula(object@splines[[1]]$call))[c(2,1,3)], collapse=""),"\n","Number of intervals:", object@number),collapse="\n" )))

####this function is used to speed up the algorithm Whenever a new number of intervals are found, this is called to
#ensure every lambda which is ahead of it in the grid is smaller than it.
arrange = function(x, y,k){
  ta=which(diff(x)>0)
  ta=ta[which(y[ta]!=0)]
  
  rms = which(ta==k)
  if(length(rms)>0){
    ta = ta[-rms]
  }
  
  if(length(ta)>0){
    if(ta==(length(x)-1)){
      x[length(x)]=0
      x[length(x)-1]=.000001
      return(x)
    }
    these = which(diff(x)>0)
    
    for(i in 1: length(these)){
      cutt = x[these][i]
      ends = which(y!=0)[which(which(y!=0)%in%which(x<cutt))][1]
      mids = these:ends
      y[mids[-c(1, length(mids))]]=0
      x[mids]=seq(from = cutt,to=x[ends], length.out = length(mids) )
      
    }
  }
  return(x)
}


###this is a function which makes sure the splines are the splines which correspond to the right value of lambda
fix.spline = function(ans, D.tid, tau, modes ){
  
  fix=1:length(ans$c.funs)
  x_loc = D.tid$x_loc
  
  for(p in 1:length(fix)){
    
    if(modes == "binomial"){
      fits=bigssg(meth~x_loc*status, type=c(x_loc="cub", status="nom"), family="binomial",
                  lambdas = ans$lambdas[fix[p]], rparm=c(x_loc=.01,status=1), data =D.tid)
    }else{
      fits=bigssa(meth~x_loc*status, type=c(x_loc="cub", status="nom"), 
                  lambdas = ans$lambdas[fix[p]], rparm=c(x_loc=.01,status=1), data =D.tid)
    } 
    
    Ds = D.tid[,c(2,3)]
    levels(Ds$status)=c(names(table(fits$xvars$status))[1],names(table(fits$xvars$status))[2])
    preds = predict(fits, Ds)
    
    if(modes=="binomial"){
      preds = preds$fitted.values
    }
    
    c.fun = preds[which(Ds$status==names(table(Ds$status))[1])]-preds[which(Ds$status==names(table(Ds$status))[2])]
    
    diffs=which(abs(c.fun)>tau)
    
    
    pos = list(-1)
    neg = list(-1)
    
    Vec = diffs
    Breaks <- c(0, which(diff(Vec) != 1), length(Vec))
    if(length(Breaks)==2){
      res = list(Vec[(Breaks[1] + 1):Breaks[2]])
    }else{
      res = sapply(seq(length(Breaks) - 1), function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
    }
    
    # # 
    # max.num = length(pos)+length(neg)
    for(i in 1:length(res)){
      if(c.fun[res[[i]][1]]<0){
        neg = c(neg, list(x_loc[res[[i]]]))
      }else{
        pos = c(pos, list(x_loc[res[[i]]]))
      }
    }
    
    pos=pos[-1]
    neg=neg[-1]
    
    lens = length(pos)+length(neg)
    
    ans$`pos ints`[fix[p]] = list(pos)
    ans$`neg ints`[fix[p]] = list(neg)
    ans$c.funs[fix[p]]=list(fits)
    ans$n.funs[fix[p]]= list(fits)
    
  }
  
  return(ans)
}


###this is the meat and potatoes of the package. It takes a list in the form output by the function prep.data and returns
###a 

###the inputs to this function are:
#D.tidy is the output of prep data

#subj is "no" or "yes" dictating whether a subject interaction is included. Be wary, true takes a great deal of memory

#tau is the value of the difference in methylation cutoff for the dmr. soon this will include a false option that uses
#standard errors instead

#high is the upper value of the smoothing spline parameter to try

#max.in requires knowledge of the algorithm to understand. After trying a value of lambda, if no new number of intervals
#is found, the value of lambda is dcremented and the process is tried again before going to the next grid values of 
#lambda. max.in manages the number of times the value of lambda is decrememnted before going to the nxt grid value of lambda

###getbw is a flag wthich dictates whether the bandwdths corresponding to each smoothing spline are calculated

###maxs also requires some knowledge of the algorithm. If all the grid values of lambda have been tried, but all of the 
#intervals have not been found, the the process is restarted but this time looking only at the intervals we havent found.
#maxs dictates the number of times the process is allowed to restart

### lam.l is a supplied grid of lambda values to try. If set to false the grid is just the an equally spaced grid 
#starting at high and going to 0

###modes is a flag that is normal or binomial and it dictates the link function

###se dictates whether the standard errors are fit


find_intervals = function(D.tidy, subj="no", tau, high=.001, max.in = 20, getbw=FALSE, maxs = 5, lam.l=FALSE, check.all = 20,modes = "binomial", se = TRUE){
  ###rad in data from prep_data and name them
  ind1 = D.tidy$x.n
  ind2 = D.tidy$x.c
  heal =D.tidy$y.n
  can = D.tidy$y.c
  
  ###since status will have two classes, we need to double the methylation and input data. 
  meth = c(heal, can)
  status = c(rep("healthy", length(heal)), rep("cancer", length(can)))
  x_loc = c(ind1,ind2)
  
  D.tid = data.frame(meth, x_loc, status)
  
  ###these loops allow for different model types. Right now the suppourted options allow for a subject variable
  ###and either normal link function (none) or a binomial link function. Subj sets whether a subject interaction is included
  ###This is where the maximum number of intervals is calculated
  if(subj!="no"){
    subj = rep(subj, 2)
    if(modes=="binomial"){
      fits=bigssg(meth~x_loc*status*subj,se.lp=se, type=c(x_loc="cub", status="nom", subj="nom"), 
                  lambdas = 0, rparm=c(x_loc=.01,status=1, subj =1), family="binomial")
      
      D.tid= as.data.frame(cbind(fits$yvar,fits$xvars$x_loc,fits$xvars$status, fits$xvars$subj))
      names(D.tid)=c("meth","x_loc", "status","subj")
      D.tid$x_loc= as.numeric(as.character(D.tid$x_loc))
      D.tid$meth=as.numeric(as.character(D.tid$meth))
      
      preds = predict(fits, D.tid[,c(2,3,4)], se.lp=se)
      preds = preds$fitted.values
    }else{
      fits=bigssa(meth~x_loc*status*subj,se.fit=se, type=c(x_loc="cub", status="nom", subj="nom"), 
                  lambdas = 0, rparm=c(x_loc=.01,status=1, subj =1))
      
      D.tidy= as.data.frame(cbind(fits$yvar,fits$xvars$x_loc,fits$xvars$status, fits$xvars$subj))
      names(D.tid)=c("meth","x_loc", "status","subj")
      D.tidy$x_loc= as.numeric(as.character(D.tid$x_loc))
      D.tidy$meth=as.numeric(as.character(D.tid$meth))
      
      preds = predict(fits, D.tid[,c(2,3,4)],se.fit=se)
      
    }
  }else{
    if(modes=="binomial"){
      fits=bigssg(meth~x_loc*status,se.lp=se, type=c(x_loc="cub", status="nom"), 
                  lambdas = 0, rparm=c(x_loc=.01,status=1), family="binomial")
      
      D.tid= as.data.frame(cbind(fits$yvar,fits$xvars$x_loc,fits$xvars$status))
      names(D.tid)=c("meth","x_loc", "status")
      D.tid$x_loc= as.numeric(as.character(D.tid$x_loc))
      D.tid$meth=as.numeric(as.character(D.tid$meth))
      
      preds = predict(fits, D.tid[,c(2,3)], se.lp = se)
      preds = preds$fitted.values
    }else{
      fits=bigssa(meth~x_loc*status,se.fit=se, type=c(x_loc="cub", status="nom"), 
                  lambdas =0 , rparm=c(x_loc=.01,status=1), data=D.tid)
      
      D.tidy= as.data.frame(cbind(fits$yvar,fits$xvars$x_loc,fits$xvars$status))
      names(D.tidy)=c("meth","x_loc", "status")
      D.tidy$x_loc= as.numeric(as.character(D.tidy$x_loc))
      D.tidy$meth=as.numeric(as.character(D.tidy$meth))
      
      preds = predict(fits, D.tidy, se.fit=se)
    }
  }
  
  
  #this is the cancer functions
  c.fun = preds[which(status==names(table(status))[1])]-preds[which(status==names(table(status))[2])]
  
  #these are the individual points which are differentially methylated
  diffs=which(abs(c.fun)>tau)
  
  #this next block of code determines the dmrs
  if(length(diffs)==0){
    print("No intervals at chosen values of tau. Lower tau and try again")
    return(0)
  }  
  
  
  
  pos = list(-1)
  neg = list(-1)
  
  Vec = diffs
  Breaks <- c(0, which(diff(Vec) != 1), length(Vec))
  if(length(Breaks)==2){
    res = list(Vec[(Breaks[1] + 1):Breaks[2]])
  }else{
    res = sapply(seq(length(Breaks) - 1), function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
  }
  # # 
  # max.num = length(pos)+length(neg)
  for(i in 1:length(res)){
    if(c.fun[res[[i]][1]]<0){
      neg = c(neg, list(res[[i]]))
    }else{
      pos = c(pos, list(res[[i]]))
    }
  }
  
  pos=pos[-1]
  neg=neg[-1]
  max.num = length(res)
  
  print(paste(c("max=",max.num), collapse=""))
  
  ###this is where the lambda grid is set
  num.int = rep(0, max.num+1)
  num.int[length(num.int)]=max.num
  len.lam = length(num.int)
  
  if(length(lam.l)==1){
    if(lam.l==FALSE){
      lam = seq(from=high,to=0, length.out=length(num.int))
      alph = mean(abs(diff(lam)))
      alph0=alph
    }
  }else{
    lam.l = lam.l
  }
  
  ##this initializes the interval lists and counters 
  int.pos.L = vector(length=len.lam, mode = "list")
  int.neg.L = int.pos.L
  
  int.pos.L[length(num.int)]= list(pos)
  int.neg.L[length(num.int)]=list(neg)
  
  can.funs = vector(length=len.lam, mode = "list")
  nor.funs = vector(length=len.lam, mode = "list")
  can.funs[length(num.int)]=list(fits)
  nor.funs[length(num.int)]= list(fits)
  
  tot.max = max.num*max.in
  j=1
  seen.one = FALSE
  countss=0
  
  ###this is the loop that contains the algorithm. 
  while((table(num.int)[1]>2)&(j<tot.max)&(countss<tot.max)){
    
    ###If we have progressed through the full lambdagrid, this value becomes positive and we enter this if loop
    if(j>length(lam)){
      nas=which(is.na(num.int))
      num.int[nas]=0
      
      j =  max(which(num.int==0))
      countss=countss+1
      
      ###if countss> maxs then the algorithm has run for too long and output whatever we have
      if(countss>maxs){
        keep = which(num.int!=0)
        if(num.int[1]==num.int[2]){
        }
        
        LL=list("tau"=tau,"pos ints"=int.pos.L[keep], "neg ints"=int.neg.L[keep], "lambdas"=lam[keep], "num.ints"= num.int[keep], "c.funs" = can.funs[keep], "n.funs"=nor.funs[keep])
        
        
        if(getbw==TRUE){
          
          k = rep(0, length(keep))
          for(i in 1:length(k)){
            k[i] =median(smooth.matrix(x.n, LL$lambdas[i], 100))
          }
          L1 = c(LL, list("bw"=k))
          return(L1)
        }else{
          return(LL)
        }
      }
      
      
      ###if not, the algorithm chooses an empty interval number and takes as the lambda guess for that values the average
      #of the two closest confirmed lambda values
      lam = arrange(lam, num.int, 1)
      lam[j]=(lam[closest.right(num.int,j-1)] + lam[(1:max.num)[closest.left(num.int,j-1)+1]])/2
      
      
      
    }
    
    
    ###this next section is the same as in the intro section but it now fits with a value of lambda from the grid 
    
    if(subj!="no"){
      subj = rep(subj, 2)
      if(modes=="binomial"){
        fits=bigssg(meth~x_loc*status*subj,se.lp=TRUE, type=c(x_loc="cub", status="nom", subj="nom"), 
                    lambdas = lam[j], rparm=c(x_loc=.01,status=1, subj =1), family="binomial")
        
        preds = predict(fits, D.tid[,c(2,3,4)])
        preds = preds$fitted.values
      }else{
        fits=bigssa(meth~x_loc*status*subj,se.fit=TRUE, type=c(x_loc="cub", status="nom", subj="nom"), 
                    lambdas = lam[j], rparm=c(x_loc=.01,status=1, subj =1))
        
        preds = predict(fits, D.tid[,c(2,3,4)])
      }
    }else{
      if(modes=="binomial"){
        fits=bigssg(meth~x_loc*status,se.lp=TRUE, type=c(x_loc="cub", status="nom"), 
                    lambdas = lam[j], rparm=c(x_loc=.01,status=1), family="binomial")
        
        preds = predict(fits, D.tid[,c(2,3)])
        preds = preds$fitted.values
      }else{
        fits=bigssa(meth~x_loc*status,se.fit=TRUE, type=c(x_loc="cub", status="nom"), 
                    lambdas = lam[j], rparm=c(x_loc=.01,status=1))
        
        preds = predict(fits, D.tid[,c(2,3)])
      }
    }
    c.fun = preds[which(status=="cancer")]-preds[which(status=="healthy")]
    
    diffs=which(abs(c.fun)>tau)
    
    if(length(diffs)==0){
      j=j+1
      next
    }
    
    pos = list(-1)
    neg = list(-1)
    
    
    Vec = diffs
    Breaks <- c(0, which(diff(Vec) != 1), length(Vec))
    if(length(Breaks)==2){
      res = list(Vec[(Breaks[1] + 1):Breaks[2]])
    }else{
      res = sapply(seq(length(Breaks) - 1), function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
    }

    for(i in 1:length(res)){
      if(c.fun[res[[i]][1]]<0){
        neg = c(neg, list(res[[i]]))
      }else{
        pos = c(pos, list(res[[i]]))
      }
    }
    
    pos=pos[-1]
    neg=neg[-1]
    
    lens = length(pos)+length(neg)
    
    
    ###check to see whether the number of intervals gotten from the grid value is new
    checkk=lens%in%num.int
    
    
    ###if so, write to list
    if((checkk==FALSE)){
      kk = length(pos)+length(neg)
      num.int[kk+1]=kk
      int.pos.L[kk+1] = list(pos)
      int.neg.L[kk+1] = list(neg)
      can.funs[kk+1]=list(fits)
      nor.funs[kk+1]= list(fits)
      lam[kk+1]=lam[j]
      lam = arrange(lam, num.int, lens)
      print(kk-1)
      j=j+1
      next
    }
    
    
    if(j!=1){
      check=(num.int[j]!=num.int[j-1]+1)&(j!=1)
      if(is.na(check)){
        j = length(num.int)+1
        next
      }
    }else{
      check = FALSE
      j=j+1
      next
    }
    
    
    ###this is the section that decrements the value lambda if the grid value didn't result in a new interval
    if(j>1){
      k=1
      
      ###choose a nearby lambda
      if((countss==0)&(j<length(lam))){
        l = abs(min(lam[j], lam[j+1]))
        
        if(is.na(l)){
          nas=which(is.na(num.int))
          num.int[nas]=0
          
          
        }
      }else{
        l = abs(lam[j])
      }
      indd= 1
      
      check.count = 1
      
      while(check){
        if(check.count>check.all){
          check = FALSE
        }
        if(indd < max.in){
          
          if(is.na(l)){
            capts = which(num.int==0)
            kk=sample(capts,1)
            l = lam[kk]
          }
          
          ###try the value of lambda
          if(subj!="no"){
            subj = rep(subj, 2)
            if(modes=="binomial"){
              fits=bigssg(meth~x_loc*status*subj,se.lp=TRUE, type=c(x_loc="cub", status="nom", subj="nom"), 
                          lambdas = l, rparm=c(x_loc=.01,status=1, subj =1), family="binomial")
              
              preds = predict(fits, D.tid[,c(2,3,4)])
              preds = preds$fitted.values
            }else{
              fits=bigssa(meth~x_loc*status*subj,se.fit=TRUE, type=c(x_loc="cub", status="nom", subj="nom"), 
                          lambdas = l, rparm=c(x_loc=.01,status=1, subj =1))
              
              preds = predict(fits, D.tid[,c(2,3,4)])
            }
          }else{
            D.tid = data.frame(meth, status, x_loc, subj)
            if(modes=="binomial"){
              fits=bigssg(meth~x_loc*status,se.lp=TRUE, type=c(x_loc="cub", status="nom"), 
                          lambdas = l, rparm=c(x_loc=.01,status=1), family="binomial")
              
              preds = predict(fits, D.tid[,c(2,3)])
              preds = preds$fitted.values
            }else{
              fits=bigssa(meth~x_loc*status,se.fit=TRUE, type=c(x_loc="cub", status="nom"), 
                          lambdas = l, rparm=c(x_loc=.01,status=1))
              
              preds = predict(fits, D.tid[,c(2,3)])
              
            }
          }
          
          c.fun = preds[which(status=="cancer")]-preds[which(status=="healthy")]
          
          diffs=which(abs(c.fun)>tau)
          
          if(length(diffs)==0){
            l=l/2
            next
          }
          pos = list(-1)
          neg = list(-1)
          
          
          Vec = diffs
          Breaks <- c(0, which(diff(Vec) != 1), length(Vec))
          if(length(Breaks)==2){
            res = list(Vec[(Breaks[1] + 1):Breaks[2]])
          }else{
            res = sapply(seq(length(Breaks) - 1), function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
          }
          

          for(i in 1:length(res)){
            if(c.fun[res[[i]][1]]<0){
              neg = c(neg, list(res[[i]]))
            }else{
              pos = c(pos, list(res[[i]]))
            }
          }
          
          
          pos=pos[-1]
          neg=neg[-1]
          temp = length(pos)+length(neg)
          print(temp)
          
          
          ###this is the section where the algorithm decides whether to decrement, increment, or keep the new 
          #interval set
          

          if(temp==(num.int[j-1]+1)){
            num.int[temp+1]=temp
            lam[temp+1]= l
            int.pos.L[temp+1]= list(pos)
            int.neg.L[temp+1]=list(neg)
            can.funs[temp+1]=list(fits)
            nor.funs[temp+1]= list(fits)
            check = FALSE
            alph=alph0
            lam = arrange(lam, num.int, temp)
            print(temp)
            check = FALSE
          }
          
          if(temp>num.int[j-1]+1){
            if(temp%in%num.int){
              
              if(temp==2){
                if(lam[temp+1]>l){
                  lam[temp+1]=l
                  lam = arrange(lam, num.int, temp)
                }
                left=lam[closest.right(num.int,j)]
                l = mean(c(l,left))
              }else{
                
                if(lam[temp+1]>l){
                  lam[temp+1]=l
                  lam = arrange(lam, num.int, temp)
                }
                
                left=lam[closest.left(num.int,j)+1]
                if(is.na(left)){
                  left=0
                }
                if(left!=0){
                  l = mean(c(l,left))
                }else{
                  l = l+alph
                }
              }
            }else{
              num.int[temp+1]=temp
              lam[temp+1]= l
              int.pos.L[temp+1]= list(pos)
              int.neg.L[temp+1]=list(neg)
              can.funs[temp+1]=list(fits)
              nor.funs[temp+1]= list(fits)
              alph = alph0
              lam = arrange(lam, num.int, temp)
              print(temp)
              check = FALSE
            }
          }else if(temp==num.int[j-1]){
            
            if(l<lam[temp+1]){
              lam[temp+1]=l
              lam = arrange(lam, num.int, temp)
            }
            
            l1 = mean(c(l,lam[closest.right(num.int,j-1)]))
            if(l1<0){
              l=l/2
            }else{
              l=l1
            }
            k=k+1
          }else if(temp<max(num.int[2:(j-1)])){
            
            if(l<lam[temp+1]){
              lam[temp+1]=l
              lam = arrange(lam, num.int, temp)
            }
            
            l1 = mean(c(l,lam[closest.right(num.int,temp)]))
            if(l1<0){
              l=l/2
            }else{
              l=l1
            }
            if((k>3)&&(k%%4==0)){
              alph = alph/2

            }
            if(l<0){
              l=abs(l)
            }
            print(temp)


            indd = indd+1
          }else{
            lam[temp+1]= l
            int.pos.L[temp+1]= "NA"
            int.neg.L[temp+1]="NA"
            check = FALSE
            can.funs[temp+1]="NA"
            nor.funs[temp+1]= "NA"
            alph = alph0
          }
        }
        len.lam = length(lam)
        check.count = check.count+1
      }
      j = j+1
    }
    
  }
  keep = which(num.int!=0)
  

  len= max(x.n)-min(x.n)
  dens = length(x.n)/len
  scaled.lam = lam/dens
  
  R1=list("tau"=tau,"pos ints"=int.pos.L[keep], "neg ints"=int.neg.L[keep], "lambdas"=lam[keep], "scaled lam"=scaled.lam[keep], "num.ints"= num.int[keep], "c.funs" = can.funs[keep], "n.funs"=nor.funs[keep])
  
  asd=fix.spline(R1, D.tid, tau, modes)
  
  LL = multiscaless(tau = R1$tau,  lambdas = R1$lambdas,Intervals=list("pos"=R1$`pos ints`,"negs"=R1$`neg ints`), splines = R1$c.funs, number=max(R1$num.ints))
  
  return(LL)
  
}
