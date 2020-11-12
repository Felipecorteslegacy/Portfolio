
library(shiny)
library(shinythemes)
library(shinydashboard)
library(zoo)
library(xts)
library(quantmod)
library(quadprog)
library(PortfolioAnalytics)
library(PerformanceAnalytics)
library(ROI)
library(ggplot2)
library(plotly)
library(dplyr)


    
    activos = c("CTXS", "QCOM", "AAPL", "IDXX", "ALGN", "AMGN", "CPRT", "FAST", 
                "CTAS", "ORLY", "EBAY", "BKNG", "NFLX", "EA", "GOOGL", "PEP", 
                "MNST", "COST", "XEL", "EXC")
    
    
    prices = read.csv("www/data_set.csv", sep = ";", dec = ",")
    
    prices = xts(prices, order.by = as.Date(rownames(prices), "%d/%m/%Y"))
    
    fechai = "2008-01-01"
    fechaf = "2020-10-30"
    periodicidad = "monthly"
    
    retornos = diff(log(prices))[-1]
    
    mu <- (colMeans(retornos))*12
    cov <- cov(retornos)*12
    var <- diag(cov)
    sigma <- sqrt(var)
    
    
    in.sample = 0.9
    
    out.sample <- 1-in.sample
    
    
    t <- length(retornos[,1])
    t.in <- round(in.sample*t, digits = 0)
    t.out <- t-t.in
    
    ret.in <- retornos[1:t.in,]
    ret.out <- retornos[t.in:t,][-1]
    
    r.sample <- list()
    r.sample[[1]] <- ret.in
    r.sample[[2]] <- ret.out
    


server <- function(input, output) {

general = reactive({        
        ## Clasificación de activos
        rf <- input$rf # rf anual
        n.act <- 20 # No. activos
        sharpe.act <- (mu-rf)/sigma
        
        ### Optimización: In.sample
        ret <<- ret.in
        short = 0   #(1: con cortos; 0: sin cortos)
        
        ## Markowitz
        
        markowitz <<- function(ret){
            ret <- ret
            mu <- mu
            cov <- cov
            n <- ncol(ret)
            T <- nrow(ret)
            # Sin restricciones en corto
            if(short == 1){
                ones <- rep(1,n)
                x <- t(mu)%*%solve(cov)%*%mu
                y <- t(mu)%*%solve(cov)%*%ones
                z <- t(ones)%*%solve(cov)%*%ones
                d <- x*z - y*y
                g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
                h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)
                rpmin <- min(mu)
                rpmax <- max(mu)
                nport <- 100
                j <- seq(rpmin,rpmax, length=nport) 
                wpo <- matrix(c(0), ncol=n, nrow=nport) 
                rpo <- matrix(c(0), nrow=nport)
                sigmapo <<- matrix(c(0), nrow=nport)
                wj <- 0
                cont <- 1
                for(i in 1:nport){
                    wj <- g + h*j[i] 
                    wpo[cont,] <- t(wj)
                    rpo[cont,] <- t(wj)%*%mu
                    sigmapo[cont,] <- sqrt(t(wj)%*%cov%*%wj)
                    cont <- cont+1
                }
                
                # PMVG
                cov_inv_1 <- solve(cov, ones) 
                wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1
                rpmvg <-mu%*%wpmvg
                sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
                
                FE <- list()
                FE[[1]] <- wpo
                FE[[2]] <- rpo
                FE[[3]] <- sigmapo
                FE[[4]] <- wpmvg
                FE[[5]] <- rpmvg
                FE[[6]] <- sigmapmvg
                return(FE)
                
            }
            # Con restricciones en corto
            else {
                # FE    
                library(quadprog)
                nport <- 100
                if(min(mu) > 0){rpmin = min(mu)+0.01}
                else{rpmin = 0.01}
                rpmax <- max(mu)-0.01
                n <- length(mu)
                j <- seq(rpmin,rpmax,length=nport)
                sigmapo <- matrix(0,nrow=nport)
                wpo <- matrix(0,nrow=nport, ncol=n)
                Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
                dvec <- rep(0,n) 
                Dmat <- 2*cov
                for(i in 1:nport){
                    bvec <- c(1,j[i],rep(0,n))
                    result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
                    wpo[i,] <- result$solution
                    sigmapo[i,] <- sqrt(result$value)
                }
                rpo <- j
                colnames(wpo) <- c(activos)
                # PMVG
                pmvg <- cbind(sigmapo,wpo)
                pmvg.sort <- pmvg[order(pmvg[,1]),]
                pmvg.sel <- cbind(pmvg.sort[1,])
                wpmvg <- cbind(pmvg.sel[2:length(pmvg.sel)])
                rownames(wpmvg) <- c(activos)
                rpmvg <- mu%*%wpmvg
                sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
                FE <- list()
                FE[[1]] <- wpo
                FE[[2]] <- rpo
                FE[[3]] <- sigmapo
                FE[[4]] <- wpmvg 
                FE[[5]] <- rpmvg 
                FE[[6]] <- sigmapmvg 
                return(FE)  
            }
        }
        
        FE <- markowitz(ret)
        wpo <<- FE[[1]]
        rpo <<- FE[[2]]
        sigmapo <<- FE[[3]]
        
        #PMVG
        wpmvg <<- FE[[4]]
        rpmvg <<- FE[[5]]
        sigmapmvg <<- FE[[6]]
        
        ## Sharpe
        
        sharpe <<- function(ret,rf){
            ret = ret
            rf = rf
            n = n.act
            sharpe.act <- (mu-rf)/sigma
            sort.act <- sharpe.act[order(-sharpe.act)]
            names <- names(sort.act)
            clasif.sharpe <- cbind(ret[,names[1:n]])
            activos <- names(clasif.sharpe)
            ret <- clasif.sharpe
            mu <- (colMeans(ret))*12
            cov <- cov(ret)*12
            var <- diag(cov)
            sigma <- sqrt(var)
            
            # Se permiten cortos
            if(short == 1){
                Er <- mu-rf 
                Z <- solve(cov,Er)  
                sumZ <- sum(Z) 
                wpt <- Z/sumZ 
                rpt <- t(wpt)%*%mu
                sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
                
                PT <- list()
                PT[[1]] <- wpt
                PT[[2]] <- rpt 
                PT[[3]] <- sigmapt
                PT[[4]] <- (rpt-rf)/sigmapt
                return(PT)
            }
            # No se permiten cortos
            else {
                library(quadprog)
                nport <- 100
                if(min(mu) > 0){rpmin = min(mu)+0.01}
                else{rpmin = 0.01}
                rpmax <- max(mu)-0.01
                
                j <- seq(rpmin,rpmax, length=nport)
                wpo <- matrix(c(0), ncol=n, nrow=nport) 
                sigmapo <- matrix(c(0), nrow=nport)
                
                Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
                dvec <- rep(0,n) 
                Dmat <- 2*cov
                for(i in 1:nport){
                    bvec <- c(1,j[i],rep(0,n))
                    result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
                    wpo[i,] <- result$solution
                    sigmapo[i,] <- sqrt(result$value)
                }
                rpo <- j
                colnames(wpo) <- c(activos)
                sharpe_port <- (rpo-rf)/sigmapo
                sharpe <- cbind(sharpe_port,wpo)
                sharpe.sort <- sharpe[order(-sharpe[,1]),]
                sharpe.sel <- cbind(sharpe.sort[1,])
                wpt <- round(cbind(sharpe.sel[2:length(sharpe.sel)]),4)
                rownames(wpt) <- c(activos)
                rpt <- mu%*%wpt
                sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
                
                PT <- list()
                PT[[1]] <- wpt
                PT[[2]] <- rpt
                PT[[3]] <- sigmapt
                PT[[4]] <- (rpt-rf)/sigmapt
                return(PT)
            }
        }
        
        PT <- sharpe(ret,rf)
        wpt <<- t(PT[[1]]) 
        rpt <<- PT[[2]] 
        sigmapt <<- PT[[3]]
        maxsharpe <<- (rpt-rf)/sigmapt
        
    
        precios <- function(activos){ 
            fechai <- fechai
            fechaf <- fechaf
            dist <- round(as.numeric((as.Date(fechaf) - as.Date(fechai))/30)-2)
            
            precios <- xts()
            for(i in 1:length(activos)){
                tmp <- Ad(getSymbols(activos[i], from=fechai, to=fechaf, 
                                     periodicity=periodicidad, auto.assign = FALSE))
                tmp <- na.approx(tmp, na.rm=FALSE)
                precios <- cbind(precios, tmp)
            }
            colnames(precios) <- activos 
            tclass(precios) <- "Date"
            if(dist < 60){
                print("Atención!! Tienes menos de 60 datos")
            }
            else {
                print("Los precios han sido importados correctamente")
            }
            return (precios)
        }
        
        
        ## Treynor
        
        
        treynor <- function(retornos,rindice,rf){
            mu <- mu
            rf <- rf
            sigma <- sigma
            nret <- length(retornos[,1])
            if(nret < 60){
                print("Error!! No hay suficientes datos")
            }
            else {
                corte <- nret-59
                ractivos60 <- cbind(retornos[corte:nret,])
                rindice60 <- cbind(rindice[corte:nret,])
                nact <- length(activos)
                betas <- matrix(0, ncol=nact)
                for(i in 1:nact){
                    lmcapm <- lm(ractivos60[,i]~rindice60)
                    coef <- lmcapm[["coefficients"]]
                    betas[i] <- coef[2]
                }
                sigmaindice <- sd(rindice)*sqrt(12)
                
                treynora <- (mu-rf)/betas
                sigma2<-sigma^2
                mu2 <- mu^2
                betas2<- betas^2
                sigmaindice2<-sigmaindice^2
                
                varerror <- abs(sigma2 - betas2*sigmaindice2)
                
                matrix <- t(rbind(treynora,t(mu),t(sigma),betas,varerror))
                sortmatrix <- matrix[order(-matrix[,1]),]
                colnames(sortmatrix) <- c("treynor","mu","sigma","betas","varerror")
                primabeta <- (sortmatrix[,"mu"]-rf)*sortmatrix[,"betas"]
                
                ratio1 <- primabeta/sortmatrix[,"varerror"]
                ratio2 <- (sortmatrix[,"betas"]^2)/sortmatrix[,"varerror"]
                
                sum1<-cumsum(ratio1)
                sum2<-cumsum(ratio2)
                
                cutoffrate<-(sigmaindice2*sum1)/(1+sigmaindice2*sum2)
                diff <- sortmatrix[,1] - cutoffrate
                cond.diff <- diff[!is.na(diff) & diff > 0]
                N.optimo <- length(cond.diff)
                cuttof.r <- cutoffrate[N.optimo]
                
                Zi <- (sortmatrix[,"betas"]/sortmatrix[,"varerror"])*(sortmatrix[,"treynor"]- cuttof.r)
                Zi <- pmax(Zi,0)
                sZi <- sum(Zi)
                witp <- Zi/sZi 
                rpot <- mu%*%witp
                sigmapot <- sqrt(t(witp)%*%cov%*%witp)
                wpot <- witp[1:N.optimo]
                
                TR <- list()
                TR[[1]] <- round(wpot,4)
                TR[[2]] <- rpot
                TR[[3]] <- sigmapot
                TR[[4]] <- sortmatrix[,"treynor"]
                TR[[5]] <- cutoffrate
                TR[[6]] <- N.optimo
                return(TR)
            }
        }
        
        indice <- c(input$Bench)
        p.indice <<- precios(indice)
        r.indice <<- diff(log(p.indice))[-1]
        
        TR <- treynor(ret,r.indice,rf)
        wpot <<- t(TR[[1]])
        
        rpot <<- TR[[2]]
        sigmapot <<- TR[[3]]
        r.treynor <- TR[[4]]
        cutoffrate <<-TR[[5]] 
        N.optimo <<- TR[[6]] 
        
        # Sortino (Semi-Varianza)
        sortino <- function(ret,h){
            # Medida de riesgo como semivarianza
            ret = ret
            rf = rf
            n.act = n.act
            semiretornos <- pmin(ret,h) 
            semicovarianzas <- cov(semiretornos)*12
            semivarianzas <- diag(semicovarianzas)
            semisigma <- sqrt(semivarianzas)
            semimu <- (colMeans(semiretornos))*12
            sortino.act <- (mu-h)/semisigma
            sort.act <- sortino.act[order(-sortino.act)]
            names <- names(sort.act)
            clasif.sortino <- cbind(retornos[,names[1:n.act]])
            
            # Partimos de la clasificación de Sortino
            s.ret <- clasif.sortino
            s.mu <- (colMeans(s.ret))*12
            semi.ret <- pmin(s.ret,h)
            semi.cov <- cov(semi.ret)*12
            
            if(short == 1){
                # Se permiten cortos
                Er <- s.mu-rf # rf: prima de los activos
                Z <- solve(semi.cov,Er)  
                sumZ <- sum(Z) 
                wps <- Z/sumZ 
                wpst <- matrix(wps,nrow=1)
                colnames(wpst) <- c(names(clasif.sortino))
                rps <- s.mu%*%t(wpst)
                sigmaps <- sqrt(wpst%*%semi.cov%*%t(wpst))
                
                SR <- list()
                SR[[1]] <- wpst
                SR[[2]] <- rps
                SR[[3]] <- sigmaps
                SR[[4]] <- (rps-rf)/sigmaps
                return(SR)
            }
            # No se permiten cortos: short=0
            else {
                ################
                library(quadprog)
                minrp=min(s.mu)+0.01
                maxrp=max(s.mu)-0.01
                nport = 100
                n <- length(s.mu)
                
                mup <- seq(minrp,maxrp,length=nport)
                sigma.po <- matrix(0,nrow=nport)
                w.po <- matrix(0,nrow=nport, ncol=n)
                
                Amat <- t(rbind(rep(1,n),s.mu,diag(1,nrow=n)))
                dvec <- rep(0,n) 
                Dmat <- 2*semi.cov
                
                for(i in 1:(nport)){
                    bvec <- c(1,mup[i],rep(0,n))
                    result <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=2)
                    w.po[i,] <- result$solution
                    sigma.po[i,] <- sqrt(result$value)
                }
                
                # Exraer portafolio tangente de Sortino
                sortino.port <- (mup-rf)/sigma.po
                m.sortino <- cbind(sortino.port,w.po)
                sortino.sort <- m.sortino[order(-m.sortino[,1]),]
                sortino.sel <- sortino.sort[1,]
                wps <- sortino.sel[2:length(sortino.sel)]
                wpst <- t(wps)
                colnames(wpst) <- c(names(clasif.sortino))
                rps <- s.mu%*%t(wpst)
                sigmaps <- sqrt(wpst%*%semi.cov%*%t(wpst))
                
                SR <- list()
                SR[[1]] <- round(wpst,4)
                SR[[2]] <- rps
                SR[[3]] <- sigmaps
                SR[[4]] <- (rps-rf)/sigmaps
                return(SR)
            }
        }
        
        h <- 0
        SR <- sortino(ret,h)
        wps <<- SR[[1]] 
        rps <<- SR[[2]] 
        sigmaps <<- SR[[3]]
        maxsortino <<- (rps-rf)/sigmaps
        
        # Omega
        
        omega <- function(ret,h){
            
            # Clasificación de los activos
            ret = ret
            n <- ncol(ret)
            omegai <- rep(NA,n)
            for(j in 1:n){
                omegai[j]= sum(pmax(ret[,j]-h,0))/sum(pmax(h-ret[,j],0))
            }
            names(omegai) <- colnames(ret)
            omega.sort <- omegai[order(-omegai)]
            names<-names(omega.sort)
            n.act
            clasif.omega <- cbind(ret[,names[1:n.act]])
            mu.o <- (colMeans(clasif.omega))*12
            cov.o <- cov(clasif.omega)*12
            var.o <- diag(cov)
            sigma.o <- sqrt(var)
            
            # Optimización
            library("ROI")
            library("ROML")
            library("ROML.portfolio")
            retornos.omega <- coredata(clasif.omega)
            if(short == 1){
                # Se permiten cortos
                m <- model()
                m$variable(portfolio, lb = -1) # the portfolio choice vector; 
                m$maximize(omega(portfolio))
                opt <- optimize(m, solver="glpk", data=list(returns = retornos.omega)) 
                pesos <- round(opt$solution[grep("portfolio", names(opt$solution))]/
                                   opt$solution[grep("z", names(opt$solution))],4)
                names(pesos) <- colnames(clasif.omega)
                wpom <- pesos
                wpomt <- matrix(wpom,nrow=1)
                colnames(wpomt) <- colnames(clasif.omega)
                rpom <- mu.o%*%t(wpomt)
                sigmapom <- sqrt(wpomt%*%cov.o%*%t(wpomt))
                ret.omega <- clasif.omega%*%t(wpomt)
                
                OM <- list()
                OM[[1]] <- wpomt
                OM[[2]] <- rpom
                OM[[3]] <- sigmapom
                OM[[4]] <- sum(pmax(ret.omega-h,0))/sum(pmax(h-ret.omega,0))
                return(OM)
            }
            else {
                # No se permiten cortos
                m <- model()
                m$variable(portfolio, lb = 0) # the portfolio choice vector; 
                m$maximize(omega(portfolio))
                opt <- optimize(m, solver="glpk", data=list(returns = retornos.omega)) 
                pesos <- round(opt$solution[grep("portfolio", names(opt$solution))]/
                                   opt$solution[grep("z", names(opt$solution))],4)
                names(pesos) <- colnames(clasif.omega)
                wpom <- pesos
                wpomt <- matrix(wpom,nrow=1)
                colnames(wpomt) <- colnames(clasif.omega)
                rpom <- mu.o%*%t(wpomt)
                sigmapom <- sqrt(wpomt%*%cov.o%*%t(wpomt))
                ret.omega <- clasif.omega%*%t(wpomt)
                
                OM <- list()
                OM[[1]] <- wpomt
                OM[[2]] <- rpom
                OM[[3]] <- sigmapom
                OM[[4]] <- sum(pmax(ret.omega-h,0))/sum(pmax(h-ret.omega,0))
                return(OM)
            }
        }
        
        h <- 0
        OM <- omega(ret,h)
        wpom <<- OM[[1]] 
        rpom <<- OM[[2]] 
        sigmapom <<- OM[[3]]
        maxomega <<- OM[[4]]
        
        
        perform.in <- function(ret){
            t <- nrow(ret)
            valor <- 100 # Valor inicial del portafolio 
            wmvg <- wpmvg
            wpt <-  wpt
            wpmvg <- t(wmvg)
            wpot <- wpot
            wps <- wps
            wpom <- wpom
            
            rport <<- matrix(0,nrow=t,ncol=6)
            colnames(rport) <<- c("PMVG","Sharpe","Treynor","Sortino",
                                 "Omega","Benchmark")
            vport <<- matrix(0,nrow=t,ncol=6)
            colnames(vport) <<- c("PMVG","Sharpe","Treynor","Sortino",
                                 "Omega","Benchmark")
            
            # Retornos
            # PMVG
            rpomvg <- ret%*%wmvg
            rport[,1] <<- rpomvg 
            
            #Sharpe
            names.sharpe <<- colnames(wpt)
            clasif.sharpe <- cbind(ret[,names.sharpe])
            rpsharpe <- clasif.sharpe%*%t(wpt)
            rport[,2] <<- rpsharpe 
            
            #Treynor
            names.treynor <<- colnames(wpot)
            clasif.treynor <- cbind(ret[,names.treynor])
            rptreynor <- clasif.treynor%*%t(wpot)
            rport[,3] <<- rptreynor
            
            # Sortino
            names.sortino <<- colnames(wps)
            clasif.sortino <- cbind(ret[,names.sortino])
            rpsortino <- clasif.sortino%*%t(wps)
            rport[,4] <<- rpsortino
            
            # Omega
            names.omega <<- colnames(wpom)
            clasif.omega  <- cbind(ret[,names.omega])
            rpomega <- clasif.omega%*%t(wpom)
            rport[,5] <<- rpomega
            
            # Benchmark
            r.benchmark <- r.indice[1:t]
            rport[,6] <<- r.benchmark
            
            # Valor del portafolio
            # PMVG
            port.mvg <<- matrix(0, nrow=t)
            port.mvg[1] <- valor
            for(i in 2:t){
                port.mvg[i] <- port.mvg[i-1]*exp(rpomvg[i-1])
            }
            vport[,1] <<- port.mvg
            
            # Sharpe
            port.sharpe <<- matrix(0, nrow=t)
            port.sharpe[1] <- valor
            for(i in 2:t){
                port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
            }
            vport[,2] <<- port.sharpe
            
            # Treynor
            port.treynor <<- matrix(0, nrow=t)
            port.treynor[1] <- valor
            for(i in 2:t){
                port.treynor[i] <- port.treynor[i-1]*exp(rptreynor[i-1])
            }
            vport[,3] <<- port.treynor
            
            # Sortino
            port.sortino <<- matrix(0, nrow=t)
            port.sortino[1] <- valor
            for(i in 2:t){
                port.sortino[i] <- port.sortino[i-1]*exp(rpsortino[i-1])
            }
            vport[,4] <<- port.sortino
            
            # Omega
            port.omega <<- matrix(0, nrow=t)
            port.omega[1] <- valor
            for(i in 2:t){
                port.omega[i] <- port.omega[i-1]*exp(rpomega[i-1])
            }
            vport[,5] <<- port.omega
            
            # Benchmark
            v.benchmark <<- matrix(0, nrow=t)
            v.benchmark[1] <- valor
            
            for(i in 2:t){
                v.benchmark[i] <- v.benchmark[i-1]*exp(r.benchmark[i-1])
            }
            vport[,6] <<- v.benchmark
            
            return(vport)
            
        }
        
        Perform.in <<- perform.in(ret.in)
        
        MDSS <<- min(rport[,"Sharpe"])
        MDSTRR <<- min(rport[,"Treynor"])
        MDSORT <<- min(rport[,"Sortino"])
        MDOME <<- min(rport[,"Omega"])
        
        perform.out <<- function(ret){
          t <- nrow(ret)
          valor <- 100 # Valor inicial del portafolio 
          wmvg <- wpmvg
          wpt <-  wpt
          wpmvg <- t(wmvg)
          wpot <- wpot
          wps <- wps
          wpom <- wpom
          
          rport <- matrix(0,nrow=t,ncol=6)
          colnames(rport) <- c("PMVG","Sharpe","Treynor","Sortino",
                               "Omega","Benchmark")
          vport <- matrix(0,nrow=t,ncol=6)
          colnames(vport) <- c("PMVG","Sharpe","Treynor","Sortino",
                               "Omega","Benchmark")
          # Retornos
          # PMVG
          rpomvg <- ret%*%wmvg
          rport[,1] <- rpomvg 
          
          #Sharpe
          names.sharpe <- colnames(wpt)
          clasif.sharpe <- cbind(ret[,names.sharpe])
          rpsharpe <- clasif.sharpe%*%t(wpt)
          rport[,2] <- rpsharpe 
          
          #Treynor
          names.treynor <- colnames(wpot)
          clasif.treynor <- cbind(ret[,names.treynor])
          rptreynor <- clasif.treynor%*%t(wpot)
          rport[,3] <- rptreynor
          
          # Sortino
          names.sortino <- colnames(wps)
          clasif.sortino <- cbind(ret[,names.sortino])
          rpsortino<- clasif.sortino%*%t(wps)
          rport[,4] <- rpsortino
          
          # Omega
          names.omega <- colnames(wpom)
          clasif.omega <- cbind(ret[,names.omega])
          rpomega <- clasif.omega%*%t(wpom)
          rport[,5] <- rpomega
          
          # Benchmark
          t.all <- length(r.indice)
          t.in <- round(in.sample*length(r.indice), digits = 0)
          t.out <- t.all-t.in
          r.benchmark <- r.indice[t.in:t.all][-1]
          rport[,6] <- r.benchmark
          
          # Valor del portafolio
          # PMVG
          port.mvg <- matrix(0, nrow=t)
          port.mvg[1] <- valor
          for(i in 2:t){
            port.mvg[i] <- port.mvg[i-1]*exp(rpomvg[i-1])
          }
          vport[,1] <- port.mvg
          
          # Sharpe
          port.sharpe <- matrix(0, nrow=t)
          port.sharpe[1] <- valor
          for(i in 2:t){
            port.sharpe[i] <- port.sharpe[i-1]*exp(rpsharpe[i-1])
          }
          vport[,2] <- port.sharpe
          
          # Treynor
          port.treynor <- matrix(0, nrow=t)
          port.treynor[1] <- valor
          for(i in 2:t){
            port.treynor[i] <- port.treynor[i-1]*exp(rptreynor[i-1])
          }
          vport[,3] <- port.treynor
          
          # Sortino
          port.sortino <- matrix(0, nrow=t)
          port.sortino[1] <- valor
          for(i in 2:t){
            port.sortino[i] <- port.sortino[i-1]*exp(rpsortino[i-1])
          }
          vport[,4] <- port.sortino
          
          # Omega
          port.omega <- matrix(0, nrow=t)
          port.omega[1] <- valor
          for(i in 2:t){
            port.omega[i] <- port.omega[i-1]*exp(rpomega[i-1])
          }
          vport[,5] <- port.omega
          
          # Benchmark
          v.benchmark <- matrix(0, nrow=t)
          v.benchmark[1] <- valor
          for(i in 2:t){
            v.benchmark[i] <- v.benchmark[i-1]*exp(r.benchmark[i-1])
          }
          vport[,6] <- v.benchmark
          
          return(vport)
        }
        
        
        Perform.out <<- perform.out(ret.out)
        
    })
    
        
  output$RPT = renderInfoBox({
        general()
        infoBox(
            "Return TP", paste(round(rpt, 3) * 100, "%"), icon = icon("chart-line"),
            color = "purple", fill = TRUE, width = 2
        )
    })
    
  output$SIGMAPT = renderInfoBox({
        general()
        infoBox(
            "Risk TP", paste(round(sigmapt, 3) * 100, "%"), icon = icon("wave-square"),
            color = "light-blue", fill = TRUE
        )
        
    })
    
  output$Sharpe = renderInfoBox({
        general()
        infoBox(
            "Sharpe TP", paste(round((rpt-input$rf)/sigmapt, 3)), icon = icon("percentage"),
            color = "yellow", fill = TRUE
        )
        
    })
    
  output$MDS = renderInfoBox({
        general()
        infoBox(
            "Max-Drawdown", paste(round(MDSS, 3)*100, "%"), icon = icon("sort-amount-down"),
            color = "fuchsia", fill = TRUE
        )
        
    })
    
  output$PIET = renderPlotly({
    general()
      wwpt = data.frame(Assets = colnames(wpt), Values = as.vector(wpt))
      wwpt = subset(wwpt, Values>0)
      
      fig <- plot_ly()
      fig <- fig %>% add_pie(data = wwpt, labels = ~Assets, values = ~Values,
                             name = "Weights TP")
      
     
      fig
      
        })
  
  output$RTR = renderInfoBox({
      general()
      infoBox(
          "Return TR", paste(round(rpot, 3) * 100, "%"), icon = icon("chart-line"),
          color = "purple", fill = TRUE, width = 2
      )
  })
  
  output$SIGMATR = renderInfoBox({
      general()
      infoBox(
          "Risk TR", paste(round(sigmapot, 3) * 100, "%"), icon = icon("wave-square"),
          color = "light-blue", fill = TRUE
      )
      
  })
  
  output$SharpeTR = renderInfoBox({
      general()
      infoBox(
          "Sharpe TR", paste(round((rpot-input$rf)/sigmapot, 3)), icon = icon("percentage"),
          color = "yellow", fill = TRUE
      )
      
  })
  
  output$MDSTR = renderInfoBox({
      general()
      infoBox(
          "Max-Drawdown", paste(round(MDSTRR, 3)*100, "%"), icon = icon("sort-amount-down"),
          color = "fuchsia", fill = TRUE
      )
      
  })
  
  output$PIETR = renderPlotly({
      wwptrey = data.frame(Assets = colnames(wpot), Values = as.vector(wpot))
      
      
      fig <- plot_ly()
      fig <- fig %>% add_pie(data = wwptrey, labels = ~Assets, values = ~Values,
                             name = "Weights TR", domain = list(row = 0, column = 0))
      
      
      fig
  })  
  
  output$RSOR = renderInfoBox({
      general()
      infoBox(
          "Return SORT", paste(round(rps, 3) * 100, "%"), icon = icon("chart-line"),
          color = "purple", fill = TRUE, width = 2
      )
  })
  
  output$SIGMASOR = renderInfoBox({
      general()
      infoBox(
          "Risk SORT", paste(round(sigmaps, 3) * 100, "%"), icon = icon("wave-square"),
          color = "light-blue", fill = TRUE
      )
      
  })
  
  output$SharpeSOR = renderInfoBox({
      general()
      infoBox(
          "Sharpe SORT", paste(round((rps-input$rf)/sigmaps, 3)), icon = icon("percentage"),
          color = "yellow", fill = TRUE
      )
      
  })
  
  output$MDSOR = renderInfoBox({
      general()
      infoBox(
          "Max-Drawdown", paste(round(MDSORT, 3)*100, "%"), icon = icon("sort-amount-down"),
          color = "fuchsia", fill = TRUE
      )
      
  })
  
  output$PIESOR = renderPlotly({
      wwpsorti = data.frame(Assets = colnames(wps), Values = as.vector(wps))
      wwpsorti = subset(wwpsorti, Values>0)
      
      
      fig <- plot_ly()
      fig <- fig %>% add_pie(data = wwpsorti, labels = ~Assets, values = ~Values,
                             name = "Weights TR", domain = list(row = 0, column = 0))
      
      
      fig
  })  
  
  output$ROME = renderInfoBox({
      general()
      infoBox(
          "Return OME", paste(round(rpom, 3) * 100, "%"), icon = icon("chart-line"),
          color = "purple", fill = TRUE, width = 2
      )
  })
  
  output$SIGMAOME = renderInfoBox({
      general()
      infoBox(
          "Risk OME", paste(round(sigmapom, 3) * 100, "%"), icon = icon("wave-square"),
          color = "light-blue", fill = TRUE
      )
      
  })
  
  output$SharpeOME = renderInfoBox({
      general()
      infoBox(
          "Sharpe OME", paste(round((rpom-input$rf)/sigmapom, 3)), icon = icon("percentage"),
          color = "yellow", fill = TRUE
      )
      
  })
  
  output$MDOMEE = renderInfoBox({
      general()
      infoBox(
          "Max-Drawdown", paste(round(MDOME, 3)*100, "%"), icon = icon("sort-amount-down"),
          color = "fuchsia", fill = TRUE
      )
      
  })
  
  output$PIEOME = renderPlotly({
      wwpome = data.frame(Assets = colnames(wpom), Values = as.vector(wpom))
      wwpome = subset(wwpome, Values>0)
      
      
      fig <- plot_ly()
      fig <- fig %>% add_pie(data = wwpome, labels = ~Assets, values = ~Values,
                             name = "Weights TR", domain = list(row = 0, column = 0))
      
      
      fig
  })  
  
  output$EFA = renderPlotly({
    general()
    
    consol = as.data.frame(cbind(sigmapo, rpo))
    index = which(consol$rpo >= as.numeric(rpmvg))
    consol = consol[index,]
    
    gg = ggplot() + geom_line(aes(x = consol$V1, y = consol$rpo)) + geom_point(aes(x = sigma, y = mu, 
                                                                                   col = activos), show.legend = F)
    gg1 = gg + theme_light() + labs(title = "Risk- Return plane 20 assets", subtitle = "20 Assets", x = "Volatility", y = "Return", caption = "Datos sacados de Yahoo Finance")
    ggplotly(gg1)
    
    
  })
  
  output$EFP = renderPlotly({
    general()
    consol = as.data.frame(cbind(sigmapo, rpo))
    index = which(consol$rpo >= as.numeric(rpmvg))
    consol = consol[index,]
    
    retportas = c(rpt, rpot, rpmvg, rpom, rps)
    sigportas = c(sigmapt, sigmapot, sigmapmvg, sigmapom, sigmaps)
    Portafolios = c("PT", "POT", "PMVG", "PO", "PS")
    
    gg0 = ggplot() + geom_line(aes(x = consol$V1, y = consol$rpo), alpha = 0.7, size = 1.5) + geom_point(aes(x = sigportas, y = retportas, col = Portafolios), size = 6) 
    
    ggt = gg0 + theme_light() + labs(title = "Risk- Return plane 5 main portfolios", subtitle = "5 main portfolios", x = "Volatility", y = "Return", caption = "Datos sacados de Yahoo Finance")
    ggplotly(ggt)
  })
  
  output$PERIN = renderPlotly({
    general()
    
    fill=c("red","blue","purple","darkgreen","darkgray","black")
    
    colors <- c("PMVG" = "blue", "Sharpe" = "red", "Treynor" = "orange", "Sortino" = "Black", "Omega" = "darkgreen", "Benchmark" = "darkgrey")
    
    Perform.in = as.data.frame(Perform.in)
    xx = as.vector(seq(1, 138, 1))
    ggd = ggplot(data = Perform.in) + geom_line(aes(x = xx, y = PMVG, color = "PMVG")) + geom_line(aes(x = xx, y = Sharpe, color = "Sharpe")) + geom_line(aes(x = xx, y = Treynor, color = "Treynor")) +geom_line(aes(x = xx, y = Sortino, color = "Sortino")) + geom_line(aes(x = xx, y = Omega, color = "Omega")) + geom_line(aes(x = xx, y = Benchmark, color = "Benchmark")) + 
      scale_color_manual(values = colors)  
    
    ggdd = ggd + theme_light() + labs(title = "Portfolios' Performance", subtitle = "In- Sample", x = "Time", y = "Valor", caption = "Datos sacados de Yahoo Finance") 
    ggplotly(ggdd) 
    
  })
  
  output$PEROUT = renderPlotly({
    general()
    
    
    fill=c("red","blue","purple","darkgreen","darkgray","black")
    
    colors <- c("PMVG" = "blue", "Sharpe" = "red", "Treynor" = "orange", "Sortino" = "Black", "Omega" = "darkgreen", "Benchmark" = "darkgrey")
    
    Perform.out = as.data.frame(Perform.out)
    xx = as.vector(seq(1, 15, 1))
    ggd = ggplot(data = Perform.out) + geom_line(aes(x = xx, y = PMVG, color = "PMVG")) + geom_line(aes(x = xx, y = Sharpe, color = "Sharpe")) + geom_line(aes(x = xx, y = Treynor, color = "Treynor")) +geom_line(aes(x = xx, y = Sortino, color = "Sortino")) + geom_line(aes(x = xx, y = Omega, color = "Omega")) + geom_line(aes(x = xx, y = Benchmark, color = "Benchmark")) + 
      scale_color_manual(values = colors)  
    
    ggdd = ggd + theme_light() + labs(title = "Portfolios' Performance", subtitle = "Out- Sample", x = "Time", y = "Valor", caption = "Datos sacados de Yahoo Finance") 
    ggplotly(ggdd) 
  
    
  })
  
  output$EM = renderPlotly({
    general()
    mu <- (colMeans(retornos))*12
    cov <- cov(retornos)*12
    var <- diag(cov)
    sigma <- sqrt(var)
    
    ## Clasificación de activos
    rf <- input$rf 
    n.act <- N.optimo 
    sharpe.act <- (mu-rf)/sigma
    
    ### Optimización
    ret <- retornos
    short <- 0   #(1: con cortos; 0: sin cortos)
    
    ## Markowitz
    FE <- markowitz(ret)
    wpo <- FE[[1]]
    rpo <- FE[[2]]
    sigmapo <- FE[[3]]
    
    #PMVG
    wmvg <- FE[[4]]
    rpmvg <- FE[[5]]
    sigmapmvg <- FE[[6]]
    
    ## Sharpe
    PT <<- sharpe(ret,rf)
    wpt <<- t(PT[[1]]) 
    rpt <- PT[[2]] 
    sigmapt <- PT[[3]]
    maxsharpe <- (rpt-rf)/sigmapt
    
    
    
    #-- Optimización inversa
    
    names <- colnames(wpt)
    retornos.clasif <<- cbind(retornos[,names[1:N.optimo]])
    
    mu.mkt <<- (colMeans(retornos.clasif)*12)-rf # Retornos de equilibrio
    cov.mkt <<- cov(retornos.clasif)*12
    var.mkt <- diag(cov.mkt)
    sigma.mkt <- sqrt(var.mkt)
    
    lambda <- (rpt-rf)/sigmapt^2
    
    w.mkt <- solve(c(lambda)*cov.mkt)%*%mu.mkt
    pi <- (cov.mkt%*%w.mkt)*c(lambda)
    
    ## Modelo Black-Litterman
    
    q <- c(0.044, -0.066, 0.121, -0.076) 
    
    P <<- rbind(c(0,1,0,0,0,0,0,0,0,0),
               c(0,0,1,0,0,0,0,0,0,0),
               c(1,0,0,0,0,0,0,0,0,0),
               c(0,0,0,0,1,0,0,0,0,0))
    
    tau <- 0.05
    k <- dim(P)[1] # nrow(P)
    Omega <- matrix(0,k,k)
    
    for(i in 1:k){
      Omega[i,i] <- (t(P[i,])%*%cov.mkt%*%P[i,])*tau
    }
    
    rp.bl <- solve(solve(tau*cov.mkt)+t(P)%*%(solve(Omega)%*%P)) %*% (solve(tau*cov.mkt)%*%pi+t(P)%*%(solve(Omega)%*%q))      
    
    w.bl <- (solve(c(lambda)*cov.mkt)%*%rp.bl)/sum((solve(c(lambda)*cov.mkt)%*%rp.bl))
    
    w.mkt = as.data.frame(w.mkt)
    colnames(w.mkt) = c("Weights")
    Names0 = rownames(w.mkt)
    w.mkt$Names = Names0
    
    
    gtm0 = ggplot(data = w.mkt) + geom_bar(aes(y = Weights, x = Names0, fill = Names0), stat = "identity")
    
    gtm0 = gtm0 + labs(title = "", x = "Assets", y = "%")
    
    ggplotly(gtm0)
    
    
    
    
    
  })
  
  output$BLM = renderPlotly({
    general()
    mu <- (colMeans(retornos))*12
    cov <- cov(retornos)*12
    var <- diag(cov)
    sigma <- sqrt(var)
    
    ## Clasificación de activos
    rf <- input$rf 
    n.act <- N.optimo 
    sharpe.act <- (mu-rf)/sigma
    
    ### Optimización
    ret <- retornos
    short <- 0   #(1: con cortos; 0: sin cortos)
    
    ## Markowitz
    FE <- markowitz(ret)
    wpo <- FE[[1]]
    rpo <- FE[[2]]
    sigmapo <- FE[[3]]
    
    #PMVG
    wmvg <- FE[[4]]
    rpmvg <- FE[[5]]
    sigmapmvg <- FE[[6]]
    
    ## Sharpe
    PT <<- sharpe(ret,rf)
    wpt <<- t(PT[[1]]) 
    rpt <- PT[[2]] 
    sigmapt <- PT[[3]]
    maxsharpe <- (rpt-rf)/sigmapt
    
    
    
    #-- Optimización inversa
    
    names <- colnames(wpt)
    retornos.clasif <<- cbind(retornos[,names[1:N.optimo]])
    
    mu.mkt <<- (colMeans(retornos.clasif)*12)-rf # Retornos de equilibrio
    cov.mkt <<- cov(retornos.clasif)*12
    var.mkt <- diag(cov.mkt)
    sigma.mkt <- sqrt(var.mkt)
    
    lambda <- (rpt-rf)/sigmapt^2
    
    w.mkt <- solve(c(lambda)*cov.mkt)%*%mu.mkt
    pi <- (cov.mkt%*%w.mkt)*c(lambda)
    
    ## Modelo Black-Litterman
    
    q <- c(0.044, -0.066, 0.121, -0.076) 
    
    P <<- rbind(c(0,1,0,0,0,0,0,0,0,0),
                c(0,0,1,0,0,0,0,0,0,0),
                c(1,0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,1,0,0,0,0,0))
    
    tau <- 0.05
    k <- dim(P)[1] # nrow(P)
    Omega <- matrix(0,k,k)
    
    for(i in 1:k){
      Omega[i,i] <- (t(P[i,])%*%cov.mkt%*%P[i,])*tau
    }
    
    rp.bl <- solve(solve(tau*cov.mkt)+t(P)%*%(solve(Omega)%*%P)) %*% (solve(tau*cov.mkt)%*%pi+t(P)%*%(solve(Omega)%*%q))      
    
    w.bl <- (solve(c(lambda)*cov.mkt)%*%rp.bl)/sum((solve(c(lambda)*cov.mkt)%*%rp.bl))
    
    
    
    w.bl = as.data.frame(w.bl)
    colnames(w.bl) = c("Weights")
    Names1 = rownames(w.bl)
    w.bl$Names = Names1
    
    gtm1 = ggplot(data = w.bl) + geom_bar(aes(y = Weights, x = Names1, fill = Names1), stat = "identity")
    
    gtm1 = gtm1 + labs(title = "", x = "Assets", y = "%")
    
    ggplotly(gtm1)
    
  })
  
  output$PerBL = renderPlotly({
    
    mu <- (colMeans(retornos))*12
    cov <- cov(retornos)*12
    var <- diag(cov)
    sigma <- sqrt(var)
    
    ## Clasificación de activos
    rf <- input$rf # rf anual
    n.act <<- N.optimo # No. activos
    sharpe.act <- (mu-rf)/sigma
    
    ### Optimización
    ret <- retornos
    short <- 0   #(1: con cortos; 0: sin cortos)
    
    ## Markowitz
    FE <- markowitz(ret)
    wpo <- FE[[1]]
    rpo <- FE[[2]]
    sigmapo <- FE[[3]]
    
    #PMVG
    wmvg <- FE[[4]]
    rpmvg <- FE[[5]]
    sigmapmvg <- FE[[6]]
    
    ## Sharpe
    
    sharpe <<- function(ret,rf){
      ret = ret
      rf = rf
      n = n.act
      sharpe.act <- (mu-rf)/sigma
      sort.act <- sharpe.act[order(-sharpe.act)]
      names <- names(sort.act)
      clasif.sharpe <- cbind(ret[,names[1:n]])
      activos <- names(clasif.sharpe)
      ret <- clasif.sharpe
      mu <- (colMeans(ret))*12
      cov <- cov(ret)*12
      var <- diag(cov)
      sigma <- sqrt(var)
      
      # Se permiten cortos
      if(short == 1){
        Er <- mu-rf 
        Z <- solve(cov,Er)  
        sumZ <- sum(Z) 
        wpt <- Z/sumZ 
        rpt <- t(wpt)%*%mu
        sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
        
        PT <- list()
        PT[[1]] <- wpt
        PT[[2]] <- rpt 
        PT[[3]] <- sigmapt
        PT[[4]] <- (rpt-rf)/sigmapt
        return(PT)
      }
      # No se permiten cortos
      else {
        library(quadprog)
        nport <- 100
        if(min(mu) > 0){rpmin = min(mu)+0.01}
        else{rpmin = 0.01}
        rpmax <- max(mu)-0.01
        
        j <- seq(rpmin,rpmax, length=nport)
        wpo <- matrix(c(0), ncol=n, nrow=nport) 
        sigmapo <- matrix(c(0), nrow=nport)
        
        Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
        dvec <- rep(0,n) 
        Dmat <- 2*cov
        for(i in 1:nport){
          bvec <- c(1,j[i],rep(0,n))
          result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
          wpo[i,] <- result$solution
          sigmapo[i,] <- sqrt(result$value)
        }
        rpo <- j
        colnames(wpo) <- c(activos)
        sharpe_port <- (rpo-rf)/sigmapo
        sharpe <- cbind(sharpe_port,wpo)
        sharpe.sort <- sharpe[order(-sharpe[,1]),]
        sharpe.sel <- cbind(sharpe.sort[1,])
        wpt <- round(cbind(sharpe.sel[2:length(sharpe.sel)]),4)
        rownames(wpt) <- c(activos)
        rpt <- mu%*%wpt
        sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
        
        PT <- list()
        PT[[1]] <- wpt
        PT[[2]] <- rpt
        PT[[3]] <- sigmapt
        PT[[4]] <- (rpt-rf)/sigmapt
        return(PT)
      }
    }
    
    
    PT <- sharpe(ret,rf)
    wpt <<- t(PT[[1]]) 
    rpt <- PT[[2]] 
    sigmapt <- PT[[3]]
    maxsharpe <- (rpt-rf)/sigmapt
    
    
    
    
    
    
    
    #-- Optimización inversa
    
    names <- colnames(wpt)
    retornos.clasif <<- cbind(retornos[,names[1:length(names)]])
    
    mu.mkt <- (colMeans(retornos.clasif)*12)-rf # Retornos de equilibrio
    cov.mkt <<- cov(retornos.clasif)*12
    var.mkt <- diag(cov.mkt)
    sigma.mkt <- sqrt(var.mkt)
    
    lambda <- (rpt-rf)/sigmapt^2
    
    w.mkt <- solve(c(lambda)*cov.mkt)%*%mu.mkt
    pi <- (cov.mkt%*%w.mkt)*c(lambda)
    
    ## Modelo Black-Litterman
    
    q <- c(0.044, -0.066, 0.121, -0.076) 
    
    P <- rbind(c(0,1,0,0,0,0,0,0,0,0),
               c(0,0,1,0,0,0,0,0,0,0),
               c(1,0,0,0,0,0,0,0,0,0),
               c(0,0,0,0,1,0,0,0,0,0))
    
    tau <- 0.05
    k <- dim(P)[1] # nrow(P)
    Omega <- matrix(0,k,k)
    
    for(i in 1:k){
      Omega[i,i] <- (t(P[i,])%*%cov.mkt%*%P[i,])*tau
    }
    
    rp.bl <- solve(solve(tau*cov.mkt)+t(P)%*%(solve(Omega)%*%P)) %*% (solve(tau*cov.mkt)%*%pi+t(P)%*%(solve(Omega)%*%q))      
    
    w.bl <- (solve(c(lambda)*cov.mkt)%*%rp.bl)/sum((solve(c(lambda)*cov.mkt)%*%rp.bl))
    
    
    
    
    w.mkt = as.data.frame(w.mkt)
    colnames(w.mkt) = c("Weights")
    Names0 = rownames(w.mkt)
    w.mkt$Names = Names0
    
    
    w.bl = as.data.frame(w.bl)
    colnames(w.bl) = c("Weights")
    Names1 = rownames(w.bl)
    w.bl$Names = Names1
    
    
    t2 <- nrow(retornos.clasif)
    valor2 <- as.numeric(input$value) # Valor inicial del portafolio 
    rport2 <- matrix(0,nrow=t2,ncol=2)
    colnames(rport2) <- c("Sharpe","Black-Litterman")
    
    
    NNN = rownames(w.bl)
    
    retornos.clasif = retornos.clasif[,NNN]
    
    w.bl$Names = NULL
    
    w.bl = as.matrix((w.bl))
    
    rportbayesiano <- retornos.clasif%*%w.bl
    
    
    
    rport2[,1]<- rportbayesiano
    
    rportsharpe<- retornos.clasif%*%t(wpt)
    rport2[,2]<- rportsharpe
    
    vport2 <- matrix(0,nrow=t2,ncol=2)
    colnames(vport2) <- c("Sharpe","Black_Litterman")
    
    port.bay <- matrix(0, nrow=t2)
    port.bay[1] <- valor2
    for(i in 2:t2){
      port.bay[i] <- port.bay[i-1]*exp(rportbayesiano[i-1])
    }
    vport2[,1] <- port.bay
    
    port.sh <- matrix(0, nrow=t2)
    port.sh[1] <- valor2
    for(i in 2:t2){
      port.sh[i] <- port.sh[i-1]*exp(rportsharpe[i-1])
    }
    vport2[,2] <- port.sh 
    
    vport2 = as.data.frame(vport2)
    
    
    colorsB <- c("Black_Litterman" = "red", "Sharpe" = "black")
    
    ggblp = ggplot(data = vport2) + geom_line(aes(x = seq(1, 153, 1) , y = Sharpe, color = "Sharpe"), size = 1.2) + geom_line(data = vport2, aes(x = seq(1, 153, 1), y = Black_Litterman, color = "Black_Litterman"), size = 1.2)
    
    ggblp = ggblp + labs(title = "", x = "Time", y = "Value") + scale_color_manual(values = colorsB) 
    
    ggplotly(ggblp)
    
    
     
})
 
  output$EFRP = renderPlotly({
    mu <- (colMeans(retornos))*12
    cov <- cov(retornos)*12
    var <- diag(cov)
    sigma <- sqrt(var)
    
    Randi = read.csv(file = "www/Random.csv")
    Randi$X = NULL
    Randi = as.matrix(Randi)
    
    RetRan = Randi%*%mu
    
    SDev = matrix(0, nrow = 1762, ncol = 1)
    
    for (i in 1:length(Randi[,1])) {
      temp = sqrt(Randi[i,]%*%cov%*%t(t(Randi[i,])))
      SDev[i,1] = temp
    }
    
    EFR = cbind(RetRan, SDev)
    
    EFR = as.data.frame(EFR)
    
    colnames(EFR) = c("Ret", "Des")
    
    Optsd = 0.1317854
    Optmean = 0.07350797
    
    gpr = ggplot(data = EFR) + geom_point(aes(x = Des, y = Ret)) + geom_point(aes(x = Optsd, y = Optmean), colour = "red")
    
    gpr = gpr + labs(title = "", x = "Risk", y = "Return")
    
    ggplotly(gpr)

    
  })
  
  output$WRP = renderPlotly({
    
    mu <- (colMeans(retornos))*12
    cov <- cov(retornos)*12
    var <- diag(cov)
    sigma <- sqrt(var)
    
    Randi = read.csv(file = "www/Random.csv")
    Randi$X = NULL
    Randi = as.matrix(Randi)
    
    RetRan = Randi%*%mu
    
    SDev = matrix(0, nrow = 1762, ncol = 1)
    
    for (i in 1:length(Randi[,1])) {
      temp = sqrt(Randi[i,]%*%cov%*%t(t(Randi[i,])))
      SDev[i,1] = temp
    }
    
    EFR = cbind(RetRan, SDev)
    
    EFR = as.data.frame(EFR)
    
    colnames(EFR) = c("Ret", "Des")
    
    Optsd = 0.1317854
    Optmean = 0.07350797
    
    
    WR = c(0.060, 0.018, 0.012, 0.000, 0.006, 0.002, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.008, 0.082, 0.002, 0.484, 0.000, 0.030, 0.060, 0.236)
    
    nomi = c("CTXS",  "QCOM",  "AAPL",  "IDXX",  "ALGN",  "AMGN",  "CPRT",  "FAST",  "CTAS",  "ORLY",  "EBAY", "BKNG",  "NFLX", "EA", "GOOGL", "PE",  "MNST",  "COST",   "XEL",   "EXC")
    
    names(WR) = nomi
    
    WR = as.data.frame(WR)
    WR$Assets = nomi
    
    WR = subset(WR, WR>0)
    
    ggrw = ggplot(data = WR) + geom_bar(aes(y = WR, x = Assets, fill = Assets), stat = "identity")
    
    ggrw = ggrw + labs(x = "Assets", y = "%")
    
    ggplotly(ggrw)
    
    
  })
  
  output$SRPA = renderInfoBox({
    general()
    
    infoBox(title = "Sharpe", value = paste(round(((0.07350797- input$rf)/0.1317854), 3), "%"), icon = icon("percentage"), fill = T, col = "yellow", width = "100%")
  })
  
    
}
