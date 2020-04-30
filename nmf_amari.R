

nmf_amari <- function(nmf_output){
  
  nmf_nrun = nmf_output$measures$nrun
  nmf_rank = nmf_output$measures$rank
  
  amari_rabd <- nmf_output$fit %>% 
    imap(.x=.,.f=function(x,name){
      dissimilarity_rab(x,name)
    })
  
  rabd_df = map_df(amari_rabd,
                   .f=function(x) x$rabd)
  
  amari_out <- list(amari_rabd = amari_rabd,
                    rabd_df = rabd_df)
  
  return(amari_out)
  
}

dissimilarity_rab <- function(x,name) {
  
  nrun = length(x)
  niters = (nrun*(nrun-1)/2)
  a = as.numeric(name)
  
  rab = data.frame(
    r = rep(a,niters),
    a =      rep(2:nrun,times=c(1:(nrun-1))),
    b =   sequence(0:(nrun-1))
  )
  
  
  dissimilarity_rab_out <- map2(.x=rab$a,.y=rab$b,
                                .f=function(a1=.x,b1=.y) {
                                  cross_correlation_d(x[[a1]],x[[b1]])
                                })
  
  rab$d = map(dissimilarity_rab_out,
              .f=function(x) return(x$d) ) %>% as.numeric() 
    
  d_rab_out = list(dissimilarity_rab_out = dissimilarity_rab_out,
                   rabd = rab)
  return(d_rab_out)
  
  
}

# # Check 2:
# a1 = rab[4,"a"]
# b1 = rab[4,"b"]
# c1 = x[[a1]]
# c2 = x[[b1]]

cross_correlation_d <- function(c1,c2){

  d1 = coef(c1) %>% t()
  d2 = coef(c2) %>% t()
  
    C <- matrix(0,nrow=ncol(d1),ncol=ncol(d2))
    
    for(j in c(1:ncol(d1))){  
      for(k in c(1:ncol(d2))) {
        C[j,k] = cor(d1[,j],d2[,k])
      }
    }
      
      sum1 = 0
      for(j in c(1:ncol(d1))){
        sum1 = (1/ncol(d1))*(1-max(C[j,],na.rm=TRUE)) + sum1
      }
      
      sum2 = 0
      for(k in c(1:ncol(d2))){
        sum2 = (1/ncol(d2))*(1-max(C[,k],na.rm=TRUE)) + sum2
      }  
      
      d = (sum1 + sum2)/2
      
      cc_d_output = list(C = C,
                         d = d)
      
      return(cc_d_output)
      
  
  
  
}



