mean_sd=function(data, x){return(sapply(x,simplify=TRUE,function(y){paste0(round(mean(data[,y],na.rm=T),2),"(",round(sd(data[,y],na.rm=T),2),")")}))}

thetaSE.eap<-function(ipar,resp.data,maxCat=5,model=1,minTheta=-4.0,maxTheta=4.0,inc=0.1,prior.dist=1,prior.mean=0.0,prior.sd=1.0,D=1.0) {
  ni<-nrow(ipar); #number of items
  nExaminees<<-nrow(resp.data); #number of examiness
  NCAT<-ipar$NCAT;
  
  theta<-seq(minTheta,maxTheta,inc);
  nq<-length(theta);
  
  if (prior.dist==1) {
    prior<-dnorm((theta-prior.mean)/prior.sd); #normal prior
  } else if (prior.dist==2) {
    prior<-exp((theta-prior.mean)/prior.sd)/(1+exp((theta-prior.mean)/prior.sd))^2; #logistic prior
  } else if (prior.dist==3) {
    prior<-rep(1,nq); #uniform prior
  }
  DISC<-ipar[["a"]];
  CB<-ipar[paste("cb",1:(maxCat-1),sep="")];
  
  prep.prob<-function(){
    pp<-array(0,c(nq,ni,maxCat));
    if (model==1) {
      for (i in 1:ni) {
        ps<-matrix(0,nq,NCAT[i]+1);
        ps[,1]<-1;
        ps[,NCAT[i]+1]<-0;
        for (k in 1:(NCAT[i]-1)) {
          ps[,k+1]<-1/(1+exp(-D*DISC[i]*(theta-CB[i,k])));
        }
        #pp[,i,1]<-1-ps[,1];
        #pp[,i,NCAT[i]]<-ps[,NCAT[i]];
        for (k in 1:NCAT[i]) {
          pp[,i,k]=ps[,k]-ps[,k+1];
        }
      }
    } else if (model==2) {
      for (i in 1:ni) {
        cb<-unlist(CB[i,]);
        cb<-c(0,cb);
        zz<-matrix(0,nq,NCAT[i]);
        sdsum<-0;
        den<-rep(0,nq);
        
        for (k in 1:NCAT[i]) {
          sdsum<-sdsum+cb[k];
          zz[,k]<-exp(D*DISC[i]*(k*theta-sdsum));
          den<-den+zz[,k];
        }
        for (k in 1:NCAT[i]) {
          pp[,i,k]<-zz[,k]/den;
        }
      }
    }
    
    return(pp);
  }
  
  pp<-prep.prob();
  
  calcEAP<-function() {
    posterior<-matrix(rep(prior,nExaminees),nExaminees,nq,byrow=T);
    
    for (i in 1:ni) {
      resp<-matrix(resp.data[,i],nExaminees,1);
      prob<-t(pp[,i,resp]);
      prob[is.na(prob)]<-1.0
      posterior<-posterior*prob;
    }
    EAP<-as.vector(posterior%*%theta/rowSums(posterior));
    SE<-as.vector(sqrt(rowSums(posterior*(matrix(theta,nExaminees,nq,byrow=T)-matrix(EAP,nExaminees,nq))^2)/rowSums(posterior)));
    return(list(theta=EAP,SE=SE))
  }
  return(data.frame(calcEAP()));
}
