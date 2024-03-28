
#' @title MARSHAL: MAize Root System Hydraulic Architecture Solver
#' @param table_data    A data frame with the CRootBox simulation results. The root system architecture should be descretize by small root segment. Needs column node1ID, node2ID, branchID, x1, y1, z1, x2, y2, z2, radius [cm], length [cm], R, G, B, time [d], type, age [d], rep.
#' @param table_cond    A data frame with the plant conductivity parameters. Needs column id, order_id [root type], order [root type name], type ("Kr" or "Kx"), x (distance from the apex)[cm], y (value for Kx [cm4 hPa-1 d-1] or Kr [cm hPa-1 d-1]).
#' @param table_soil    A data frame with the soil humidity profile. Needs column: id, z [cm], psi [hPa].
#' @param hetero    Do we need to compute the uptake in heterogeneous soil? Default: hetero = TRUE
#' @param Psi_collar Water potential at the collar (boundary condition). Default: -15.000 [hPa]
#' @param verbatim If true, outputs the different step of the procedure
#' @param soil_param A data frame with the soil hydraulic parameter data(n,alpha,Ksat,lambda,Q_r,Q_s)
#' @keywords root, water
#'
#' @description Get the hydraulic properties of a given root system. Needs the root architecture from CRootBox, hydraulic properties and soil conditions.
#' @author Felicien Meunier, Adrien Heymans

require(dplyr)
getSUF <- function(table_data = NULL,
                   table_cond = NULL,
                   table_soil = NULL,
                   hetero = TRUE,
                   Psi_collar = -15000,
                   soil_param = NULL,
                   verbatim = FALSE){
  ####################################################
  #	Calculates Couvreur Macroscopic parameters   #
  ####################################################
  # Felicien Meunier, 06/2017
  # Adrien Heymans, 08/2019 - 08/2023

  if(verbatim) message("Preprocess data")

  if(length(table_data$K_type[1]) > 0){
    # If root type has been modified as a function of the radius
    radtype = T
    }else{
     radtype = F
  }

  ####################################################
  # Input data
  ####################################################
  prev <- table_data$node1ID 	  # mother segment
  l <- table_data$length    	  # segment length

  if(length(l[l == 0]) > 0){
    #message(paste0(length(l[l == 0])," root segment are 0 mm long"))
  }
  l[l == 0] <- 10e-9
  r <- table_data$radius 	  # segment radius in cm

  if(length(r[r== 0]) > 0){
   # message(paste0(length(r[r == 0])," root segment have a radius of 0 mm "))
  }
  r[r == 0] <- 10e-9


  z2 <- table_data$z2           # z-position
  z = z2

  if(radtype){
    order <- table_data$K_type
  }else{
    order <- table_data$type
  }

  seg_age <- table_data$time     # segment age
  seg_age <- max(table_data$time) - seg_age

  Nseg=length(l)      	  # Total number of segment

  Psi_sr_homogeneous <- -100      # Homogeneous soil-root potential
  Psi_sr_heterogeneous <- -3000   # Heterogeneous soil-root potential

# In addition (sept 2020), input for Ki,eq
# equation from: van Lier, 2006; Meunier, 2018, van Genuchten, 1980
  if(verbatim) message("Addition of rhizosphere")
  RLD <- table_data%>%
  mutate(z_reso = round(z2/2)*2)%>%
    dplyr::group_by(z_reso)%>%
    dplyr::summarise(rld = sum(length)/(2*1250))%>% # 8 plants/10 000 cm2 and 2 cm depth
    ungroup() %>%
    mutate(rm = sqrt(1/(pi*rld))) # radius bulk to inside of the root

  # table_data$rbulk <- 1000
  table_data$z_reso <- round(table_data$z2/2)*2
  table_data <- left_join(table_data, RLD%>%select("rm", "z_reso"), by= "z_reso")%>%select(-"z_reso")


  table_data$rho <-  ifelse(table_data$rm > r, table_data$rm/r, 1.00001)
  if(length(which(table_data$rm < r)) > 0){
    print("r_bulk < root radius")
  }
  # Meunier et al. 2018 Measuring and Modeling Hydraulic Lift of Lolium multiflorum
  # Using Stable Water Isotopes
  table_data$Be <- (2*(1-table_data$rho^2))/((-2*table_data$rho^2)*(log(table_data$rho, base = exp(1))-(1/2))-1) # dimenssion less geometric param
  # print("geometric param on root")
  if(is.null(soil_param)){
    soil_param <- data.frame (n = 1.89, alpha = 0.075,  Ksat = 25, lambda = 0.5,Q_r = 0.078, Q_s = 0.43)
  }
  soil_param$m <- 1 - 1/soil_param$n # depending on the hypothesis it could be also "m = 1 - 2/n"

  table_soil$theta_h = soil_param$Q_r+(soil_param$Q_s-soil_param$Q_r)/(1+abs(soil_param$alpha*table_soil$psi)^soil_param$n)^soil_param$m
  table_soil$theta_h[table_soil$theta_h >= soil_param$Q_s] = soil_param$Q_s
  table_soil$Se = (table_soil$theta_h-soil_param$Q_r)/(soil_param$Q_s -soil_param$Q_r)
  # table_soil$Se <- 1/(1+(abs(soil_param$alpha*table_soil$psi))^soil_param$n)^soil_param$m #

  # Van Genuchten 1980
  table_soil$Ksoil <- soil_param$Ksat*table_soil$Se^soil_param$lambda*(1-(1-table_soil$Se^(1/soil_param$m))^soil_param$m)^2
  table_soil$z_reso <- table_soil$z
  table_data$z_reso <- round(table_data$z2/2)*2
  table_data <- left_join(table_data, table_soil%>%select("psi","Se", "Ksoil", "z_reso"), by= "z_reso")
  ####################################################
  # Interpolates kr,kx functions
  if(verbatim) message("Interpolate Kr Kx function")
  order_uni=unique(order)

  kr=matrix(0,Nseg,1) # radial conductivity of the segments
  kx=matrix(0,Nseg,1) # Axial conductance of the segments

  # Linear interpolation
  for ( i in 1:length(order_uni)) {

    pos = is.element(order,order_uni[i])
    od <- order_uni[i]

    x = table_cond$x[table_cond$order_id == od & table_cond$type == "kr"]
    y = table_cond$y[table_cond$order_id == od & table_cond$type == "kr"]
    x <- c(x, 5000)
    y <- c(y, y[length(y)])

    xout = data.frame(seg_age[pos])
    temp=data.frame(approx(x,y,xout[,1]))
    kr[pos]=temp[,2]

    x = table_cond$x[table_cond$order_id == od & table_cond$type == "kx"]
    y = table_cond$y[table_cond$order_id == od & table_cond$type == "kx"]
    x <- c(x, 5000)
    y <- c(y, y[length(y)])
    temp=data.frame(approx(x,y,xout[,1]))
    kx[pos]=temp[,2]
  }

  table_data$kr <- kr[,1]
  table_data$kx <- kx[,1]

  table_data$Kieq <- (table_data$kr * table_data$Ksoil *  2*pi* table_data$length* r * table_data$Be)/( table_data$Be* table_data$Ksoil + table_data$radius * table_data$kr)
  kr_eq <- matrix(0,Nseg,1)
  kr_eq[,1] <- table_data$Kieq/(2*pi*table_data$length*r)
  kr_eq[is.na(kr_eq[,1]),1] <- 1E-15
  # Combination of hydraulics and geomitric properties
  kappa=sqrt(2*pi*r*kr*kx)  # kappa
  kappa_eq = sqrt(2*pi*r*kr_eq*kx)
  tau=sqrt(2*pi*r*kr/kx)    # tau
  tau_eq = sqrt(2*pi*r*kr_eq/kx)

  mse_kr = (1/Nseg)*sum((kr_eq-kr)^2)
  # print(paste0("MSE kr kr_eq = ",mse_kr))

  # -----------------------
  # -----------------------
  # -----------------------
  #### HOMOGENEOUS CONDITIONS
  # -----------------------
  # -----------------------
  # -----------------------

  if(verbatim) message(">> Homogeneous conditions")
  Psi_sr= Psi_sr_homogeneous * matrix(1,Nseg,1) # Soil-root potential for each segment

  ####################################################
  # Build Matrices
  if(verbatim) message("Build matrix A")
  A = Matrix(c(0),nrow=Nseg+1,ncol=Nseg+1,sparse = TRUE) # Matrix A sparse
  A_eq = Matrix(c(0),nrow=Nseg+1,ncol=Nseg+1,sparse = TRUE)

  j <- 1:Nseg
  i <- prev

  rows <- i+1
  columns <- i+1
  values=-kappa/sinh(tau*l)-kappa*tanh(tau*l/2)
  values_eq=-kappa_eq/sinh(tau_eq*l)-kappa_eq*tanh(tau_eq*l/2)

  rows=c(rows,j+1)
  columns=c(columns,i+1)
  values=c(values,kappa/sinh(tau*l))
  values_eq=c(values_eq,kappa_eq/sinh(tau_eq*l))

  rows=c(rows,i+1)
  columns=c(columns,j+1)
  values=c(values,-kappa*tanh(tau*l/2)+kappa/tanh(tau*l))
  values_eq=c(values_eq,-kappa_eq*tanh(tau_eq*l/2)+kappa_eq/tanh(tau_eq*l))

  rows=c(rows,j+1)
  columns=c(columns,j+1)
  values=c(values,-kappa/tanh(tau*l))
  values_eq=c(values_eq,-kappa_eq/tanh(tau_eq*l))
  x=mapply(values,FUN=as.numeric)
  x_eq=mapply(values_eq,FUN=as.numeric)

  A <- sparseMatrix(rows, columns, x = x) # Assignates values to specific locations
  A_eq <- sparseMatrix(rows, columns, x = x_eq)
  a <- A[-1,-1]				    # a matrix = A without the first line and column
  a_eq <- A_eq[-1,-1]

  # Build Matrix B
  if(verbatim) message("Build matrix B")
  B <- Matrix(c(0),nrow=Nseg+1,ncol=1,sparse = TRUE) # Matrix B sparse
  B_eq <- Matrix(c(0),nrow=Nseg+1,ncol=1,sparse = TRUE)

  rows <- i+1;
  columns <- matrix(1,Nseg,1)
  values <- -Psi_sr*kappa*tanh(tau*l/2)
  values_eq <- -Psi_sr*kappa_eq*tanh(tau_eq*l/2)

  rows <- c(rows,j+1)
  columns <- c(columns,matrix(1,Nseg,1))
  values <- c(values,-Psi_sr*kappa*tanh(tau*l/2))
  values_eq <- c(values_eq,-Psi_sr*kappa_eq*tanh(tau_eq*l/2))

  x <- mapply(values,FUN=as.numeric)
  x_eq=mapply(values_eq,FUN=as.numeric)
  B <- sparseMatrix(rows, columns, x = x) # Assignates values to specific locations
  B_eq <- sparseMatrix(rows, columns, x = x_eq)

  b <- B[-1] # b matrix = B without the first line
  b_eq <- B_eq[-1]

  prev_collar <- c(prev==0)

  b[prev_collar] <- b[prev_collar] - (Psi_collar * (kappa[prev_collar] / sinh(tau[prev_collar] * l[prev_collar])))
  b_eq[prev_collar] <- b_eq[prev_collar] - (Psi_collar * (kappa_eq[prev_collar] / sinh(tau_eq[prev_collar] * l[prev_collar])))

  ####################################################
  # Compute solution
  if(verbatim) message("Compute solution")
  X <- solve(a,b) 		 # a\b
  # print("solved a*b")
  X_eq <- solve(a_eq,b_eq)
  # print("solved a_eq*b_eq")

  Psi_basal <- X  		 # Solution = Psi_basal
  Psi_basal_eq <- X_eq
  prev_temp <- prev
  prev_temp[prev==0] <- 1;
  Psi_proximal <- Psi_basal[prev_temp] # Psi_proximal = Psi_basal of the mother segment
  Psi_proximal_eq <- Psi_basal_eq[prev_temp]
  Psi_proximal[prev_collar] <- Psi_collar;
  Psi_proximal_eq[prev_collar] <- Psi_collar
  Jr <- 2*kappa*tanh(tau*l/2)*(Psi_sr-(Psi_proximal+Psi_basal)/2) # Total radial flow
  Jr_eq <- 2*kappa_eq*tanh(tau_eq*l/2)*(Psi_sr-(Psi_proximal_eq+Psi_basal_eq)/2) # Total radial flow
  Jxl <- kappa*((Psi_basal-Psi_sr)/sinh(tau*l)-(Psi_proximal-Psi_sr)/tanh(tau*l)) # Axial flow at the top of the segments
  Jxl_eq <- kappa_eq*((Psi_basal_eq-Psi_sr)/sinh(tau_eq*l)-(Psi_proximal_eq-Psi_sr)/tanh(tau_eq*l))

  remove(a, b, A, B, a_eq, b_eq, A_eq, B_eq)

  # Macroscopic solution
  if(verbatim) message("Macroscopic solution")

  Tpot=sum(Jr) 		 # Actual transpiration
  Tpot_eq=sum(Jr_eq) 		 # Actual transpiration
  SUF=Jr/Tpot  		 # SUF = normalized uptake
  SUF_eq=Jr_eq/Tpot_eq
  Krs=Tpot/abs(Psi_sr_homogeneous-Psi_collar) # Total root system conductance
  Ksrs=Tpot_eq/abs(Psi_sr_homogeneous-Psi_collar)

  SUF[SUF < 0] <- 10e-10
  SUF_eq[SUF_eq < 0] <- 10e-10

  Tact <- Tpot
  Tact_eq <- Tpot_eq

  # -----------------------
  # -----------------------
  # -----------------------
  #### HETEROGENOUS CONDITIONS
  # -----------------------
  # -----------------------
  # -----------------------


  if(hetero){
    if(verbatim) message(">> Heterogeneous conditions")

    #table_soil <- rbind(data.table(id = 0, z = 100, psi = table_soil$psi[1]), table_soil)
    #table_soil <- rbind(table_soil, data.table(id = nrow(table_soil)+1, z = -1000, psi = table_soil$psi[nrow(table_soil)]))

    Psi_sr <- data.frame(approx(table_soil$z, table_soil$psi, z))[,2]

    ####################################################
    # Build Matrices
    A = Matrix(c(0),nrow=Nseg+1,ncol=Nseg+1,sparse = TRUE) # Matrix A sparse
    A_eq = Matrix(c(0),nrow=Nseg+1,ncol=Nseg+1,sparse = TRUE)

    j <- 1:Nseg
    i <- prev

    rows <- i+1
    columns <- i+1
    values=-kappa/sinh(tau*l)-kappa*tanh(tau*l/2)
    values_eq = -kappa_eq/sinh(tau_eq*l)-kappa_eq*tanh(tau_eq*l/2)

    rows=c(rows,j+1)
    columns=c(columns,i+1)
    values=c(values,kappa/sinh(tau*l))
    values_eq=c(values_eq,kappa_eq/sinh(tau_eq*l))

    rows=c(rows,i+1)
    columns=c(columns,j+1)
    values=c(values,-kappa*tanh(tau*l/2)+kappa/tanh(tau*l))
    values_eq=c(values_eq,-kappa_eq*tanh(tau_eq*l/2)+kappa_eq/tanh(tau_eq*l))

    rows=c(rows,j+1)
    columns=c(columns,j+1)
    values=c(values,-kappa/tanh(tau*l))
    values_eq=c(values_eq,-kappa_eq/tanh(tau_eq*l))
    x=mapply(values,FUN=as.numeric)
    x_eq=mapply(values_eq,FUN=as.numeric)

    if(verbatim) message("Build matrix A")
    A <- sparseMatrix(rows, columns, x = x) # Assignates values to specific locations
    a = A[-1,-1]				    # a matrix = A without the first line and column

    A_eq <- sparseMatrix(rows, columns, x = x_eq)
    a_eq = A_eq[-1,-1]

    # Build Matrix B
    if(verbatim) message("Build matrix B")
    B = Matrix(c(0),nrow=Nseg+1,ncol=1,sparse = TRUE) # Matrix B sparse
    B_eq = Matrix(c(0),nrow=Nseg+1,ncol=1,sparse = TRUE)

    rows=i+1;
    columns=matrix(1,Nseg,1)
    values=-Psi_sr*kappa*tanh(tau*l/2)
    values_eq=-Psi_sr*kappa_eq*tanh(tau_eq*l/2)

    rows=c(rows,j+1)
    columns=c(columns,matrix(1,Nseg,1))
    values=c(values,-Psi_sr*kappa*tanh(tau*l/2))
    values_eq=c(values_eq,-Psi_sr*kappa_eq*tanh(tau_eq*l/2))

    x = mapply(values,FUN=as.numeric)


    B <- sparseMatrix(rows, columns, x = x) # Assignates values to specific locations

    b=B[-1] # b matrix = B without the first line

    x_eq = mapply(values_eq,FUN=as.numeric)
    B_eq <- sparseMatrix(rows, columns, x = x_eq)
    b_eq=B_eq[-1]

    prev_collar=(prev==0)

    b[prev_collar] = b[prev_collar] - (Psi_collar * (kappa[prev_collar] / sinh(tau[prev_collar] * l[prev_collar])))
    b_eq[prev_collar] = b_eq[prev_collar] - (Psi_collar * (kappa_eq[prev_collar] / sinh(tau_eq[prev_collar] * l[prev_collar])))

    ####################################################
    # Compute solution
    if(verbatim) message("Compute solution")
    X=solve(a,b) 		 # a\b
    X_eq=solve(a_eq,b_eq)
    Psi_basal=X  		 # Solution = Psi_basal
    Psi_basal_eq=X_eq
    prev_temp=prev
    prev_temp[prev==0]=1;
    Psi_proximal=Psi_basal[prev_temp] # Psi_proximal = Psi_basal of the mother segment
    Psi_proximal_eq=Psi_basal_eq[prev_temp]
    Psi_proximal[prev_collar]=Psi_collar
    Psi_proximal_eq[prev_collar]=Psi_collar

    Jr = 2*kappa*tanh(tau*l/2)*(Psi_sr-(Psi_proximal+Psi_basal)/2) # Total radial flow
    Jr_eq = 2*kappa_eq*tanh(tau_eq*l/2)*(Psi_sr-(Psi_proximal_eq+Psi_basal_eq)/2) # Total radial flow
    Jxl = kappa*((Psi_basal-Psi_sr)/sinh(tau*l)-(Psi_proximal-Psi_sr)/tanh(tau*l)) # Axial flow at the top of the segments
    Jxl_eq = kappa_eq*((Psi_basal_eq-Psi_sr)/sinh(tau_eq*l)-(Psi_proximal_eq-Psi_sr)/tanh(tau_eq*l))

    Tact=sum(Jr) 		 # Actual transpiration
    Tact_eq=sum(Jr_eq)

    remove(a, b, A, B, a_eq, b_eq, A_eq, B_eq)

  }


  ####################################################
  if(verbatim) message("Export data")
  return(list(suf=log10(SUF),
              suf1 = SUF,
              suf_eq = SUF_eq,
              kr = log10(kr),
              kr_eq = log10(kr_eq),
              kx = log10(kx),
              tact = Tact,
              tact_eq = Tact_eq,
              tpot = Tpot,
              tpot_eq = Tpot_eq,
              krs = Krs,
              ksrs = Ksrs,
              jr = Jr,
              jr_eq = Jr_eq,
              psi = Psi_basal,
              psi_eq = Psi_basal_eq,
              jxl = Jxl,
              jxl_eq = Jxl_eq,
              psi_soil = Psi_sr,
              psi_collar = Psi_collar))
}
