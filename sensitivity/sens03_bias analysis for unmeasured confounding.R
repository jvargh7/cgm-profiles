e_value_ci = function(rr_beta,rr_lci,rr_uci){
  
  
  if(rr_beta<1){
    rr_est = 1/rr_beta
    e_val_est = rr_est + sqrt(rr_est*(rr_est - 1))
    rr_uci_est = ifelse(rr_uci >=1, rr_uci,1/rr_uci)
    
    e_val_ci = ifelse(rr_uci >= 1, 1, rr_uci_est + sqrt(rr_uci_est*(rr_uci_est-1)))
    
    
  }
  
  if(rr_beta > 1){
    rr_est = rr_beta
    e_val_est = rr_est + sqrt(rr_est*(rr_est - 1))
    
    e_val_ci = ifelse(rr_lci <= 1,1,rr_lci + sqrt(rr_lci*(rr_lci-1)))
  }
  
  return(paste0(round(e_val_est,2),"; CI: ",round(e_val_ci,2)))
  
}


e_value_ci(1.44,1.20,1.73)
e_value_ci(1.33,1.11,1.58)


e_value_ci(2.84,1.65,4.88)
e_value_ci(1.39,0.78,2.45)

e_value_ci(1.49,1.24, 1.79)
e_value_ci(1.33,1.12, 1.58)


e_value_ci(1.29,1.11,1.49)
e_value_ci(1.31,1.14,1.51)

e_value_ci(1.58, 1.25, 1.98)
e_value_ci(1.37, 1.10, 1.71)

e_value_ci(1.65,1.18,2.31)
e_value_ci(1.88,1.37,2.58)
