V 0.2.0 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10890110.svg)](https://doi.org/10.5281/zenodo.10890110)

V 0.1.0 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2391555.svg)](https://doi.org/10.5281/zenodo.2391555)

# marshal
R package for MARSHAL
The MAIZE ROOT SYSTEM HYDRAULIC ARCHITECTURE SOLVER

## Description
MARSHAL is a maize root system hydraulic architecture solver that combines the root architecture model CRootBox (Schnepf et al. 2018) with the method for solving water flow in RSHA of Meunier et al. (2017) and with the method for computing macroscopic parameter of Couvreur et al. (2012).

The model computes water flow at the level of root segments, quantifies the contribution of the water flows of each of the root segments, and predicts the whole conductivity of the root system (Krs) [cm3.hPa-1.d-1] and the potential transpiration, as well as the actual one [cm3 d-1].

## Install
To load the package, you can use the command:

```{r}
    install.packages("devtools")
    library(devtools)
    install_github("MARSHAL-ROOT/marshal")
    library(marshal)
```

## How to use it

```{r}
# Load conductivity functions (kr = f(age), Kx = f(age)) 
conductivities <- read_excel("www/conductivities.xlsx")

# Load root architecture from CPlantBox
root <- data.table::fread("./rootsystem.txt", header = TRUE)

# Create soil water profile
soil = data.frame(id=1:101,
                  z = sort(seq(-100,0,1), decreasing = TRUE),
                  value = 1,
                  psi = sort(seq(-400,-300), decreasing = TRUE))

# Set the soil hydraulic function parameters
sandy_loam = data.frame(n = 1.89, alpha = 0.075,  Ksat = 25, lambda = 0.5,Q_r = 0.078, Q_s = 0.43)

###############
# Run MARSHAL #
###############

macro_hydro = getSUF(table_data = root, 
                     table_cond = conductivities,
                     table_soil = soil, 
                     hetero = TRUE,
                     Psi_collar = -15000,
                     soil_param = sandy_loam)

# Merge output of MARSHAL on specific root segment

root$suf <- as.vector(hydraulics$suf)
root$jr <- as.vector(hydraulics$jr)
root$jr_eq <- as.vector(hydraulics$jr_eq)
root$psi <- as.vector(hydraulics$psi)
root$jxl <- as.vector(hydraulics$jxl)

# Get Whole root system conductance
Krs <-  hydraulics$krs # cm4 hPa-1 d-1
  
RLDWU <- root%>% # gather information by layer
      mutate(rz2 = round((z1+z2)/2))%>%
      dplyr::group_by(rz2)%>%
      dplyr::summarise(su = sum(suf), jr = sum(jr), jx = sum(jxl))%>%
      ungroup()

# Actual root water uptake [cm3 d-1]
Q_dou <- rev(RLDWU$jr)

# Standard uptake fraction [-]
SUF = rev(RLDWU$su)
Tact = sum(Q_dou)

# Soil-root water potential [hPa]
Hsr <- soil$psi[soil$z >= min(RLDWU$rz2) & soil$z <= max(RLDWU$rz2)] # 

# Couvreur et al. 2012 --> Kcomp
if(length(unique(Hsr))> 1){
    up = tibble(s = SUF,j = Q_dou,h = Hsr)%>%
            mutate(Hseqi = s*h, as = s*Tact_eq, upper = j - as)
    Hseq <- sum(up$Hseqi,na.omit = TRUE)
    up <- up%>%
            mutate(bel = h -Hseq, belo = (bel*s)^(-1))
    Kcomp = up$upper%*%(up$belo)
}else{
    Hseq <- Hsr[1]
    Kcomp = Krs}

```

## How to Cite

Félicien Meunier, Adrien Heymans, Xavier Draye, Valentin Couvreur, Mathieu Javaux, Guillaume Lobet, MARSHAL, a novel tool for virtual phenotyping of maize root system hydraulic architectures, in silico Plants, Volume 2, Issue 1, 2020, diz012, https://doi.org/10.1093/insilicoplants/diz012


### link to the article
[MARSHAL, a novel tool for virtual phenotyping of maize root system hydraulic architectures
Félicien Meunier, Adrien Heymans, Xavier Draye, Valentin Couvreur, Mathieu Javaux, Guillaume Lobet
in silico Plants, 2019](https://doi.org/10.1093/insilicoplants/diz012)

## Licence

MARSHAL is released under a Apache licence 2.0.
