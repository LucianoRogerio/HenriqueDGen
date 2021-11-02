
Shannon.Weaver.QT <- function(dados, min = NULL, max = NULL, formula = 1){
if(formula == 1){
dados <- as.numeric(dados)
P1 <- (sum(dados <= (mean(dados, na.rm = T) - sd(dados, na.rm = T)), na.rm = T))/(
    length(dados) - sum(is.na(dados)))
P2 <- (sum(dados <= (mean(dados, na.rm = T) - (sd(dados, na.rm = T)/2)), na.rm = T) 
   - sum(dados < (mean(dados, na.rm = T) - sd(dados, na.rm = T)), na.rm = T))/(
    length(dados) - sum(is.na(dados)))
P3 <- (sum(dados <= (mean(dados, na.rm = T)), na.rm = T) - sum(dados < (mean(dados, na.rm = T) 
   - sd(dados, na.rm = T)/2), na.rm = T))/(length(dados) - sum(is.na(dados)))
P4 <- (sum(dados <= (mean(dados, na.rm = T) + sd(dados, na.rm = T)/2), na.rm = T) 
   - sum(dados < (mean(dados, na.rm = T)), na.rm =T))/(length(dados) - sum(is.na(dados)))
P5 <- (sum(dados <= (mean(dados, na.rm = T) + sd(dados, na.rm = T)), na.rm = T) 
   - sum(dados < (mean(dados, na.rm = T) + sd(dados, na.rm = T)/2), na.rm = T))/(
    length(dados) - sum(is.na(dados)))
P6 <- (length(dados) - sum(dados < (mean(dados, na.rm = T) + sd(dados, na.rm = T)), na.rm = T)
   - sum(is.na(dados)))/(
    length(dados) - sum(is.na(dados)))

Ps <- rbind(P1, P2, P3, P4, P5, P6)
Ps[Ps == 0] <- (10^-50)
Agalinha <- -sum(Ps*log(Ps))/log(6)
return(Agalinha)} else if(formula == 2){
  dados <- as.numeric(dados)
  
  minAll <- min(dados, na.rm = T)
  maxAll <- max(dados, na.rm = T)
  ClasApt <- (maxAll - minAll)/6
  P1 <- (sum(dados <= (minAll + 1*ClasApt), na.rm = T))/(length(dados) - sum(is.na(dados)))
  
  P2 <- (sum(dados <= (minAll + 2*ClasApt), na.rm = T)
         - sum(dados < (minAll + 1*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))
  
  P3 <- (sum(dados <= (minAll + 3*ClasApt), na.rm = T)
         - sum(dados < (minAll + 2*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))
  
  P4 <- (sum(dados <= (minAll + 4*ClasApt), na.rm = T)
         - sum(dados < (minAll + 3*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))
  
  P5 <- (sum(dados <= (minAll + 5*ClasApt), na.rm = T)
         - sum(dados < (minAll + 4*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))
  
  P6 <- (sum(dados <= (maxAll), na.rm = T)
         - sum(dados < (minAll + 5*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))
  
  Ps <- rbind(P1, P2, P3, P4, P5, P6)
  Ps[Ps == 0] <- (10^-50)
  Agalinha <- -sum(Ps*log(Ps))/log(6)
  return(Agalinha)
} else{cat("This fuction offers two options of class delimitaion
limits to estimate Shannon-Weaver index.
Please choose '1' to use mean and sd to delimited Classes or
choose '2' to use the amplitude to delimited the classes")}}

Shannon.Weaver.QL <- function(dados, nclass){
  pr <- table(dados)
  pr <- pr[!is.na(pr)]
  pr2 <- pr/sum(pr)
  Agalinha <- -sum(pr2*log(pr2))/log(nclass)
return(Agalinha)}

### Option for Quantitative

dados <- DadosHen %>% dplyr::select(AntracnoseMax, ManchaPardaMax, ManchaBrancaMax, QueimaFolhasMax)
