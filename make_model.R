make_model_yk <-
  function( TmbData,
            Version,
            RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0),
            Method="Mesh",
            Npool=0,
            ConvergeTol=1,
            Use_REML=FALSE,
            loc_x=NULL,
            Parameters="generate",
            Random="generate",
            Map="generate",
            DiagnosticDir=NULL,
            TmbDir=system.file("executables",package="VAST"),
            RunDir=getwd(),
            CompileDir=TmbDir,
            build_model=TRUE ){
    #=============#
    Npool=0
    ConvergeTol=1
    Use_REML=FALSE
    loc_x=NULL
    Parameters="generate"
    Random="generate"
    Map="generate"
    DiagnosticDir=NULL
    TmbDir=system.file("executables",package="VAST")
    RunDir=getwd()
    CompileDir=TmbDir
    build_model=TRUE 
    #=============#
    TmbData = TmbData
    RunDir = DateFile
    Version = Version
    RhoConfig = RhoConfig
    loc_x = Spatial_List$loc_x
    Method = Spatial_List$Method
    #=============#
    
    # Extract Options and Options_vec (depends upon version)
    if( all(c("Options","Options_vec") %in% names(TmbData)) ){
      Options_vec = TmbData$Options_vec
      Options = TmbData$Options
    }
    if( "Options_list" %in% names(TmbData) ){
      Options_vec = TmbData$Options_list$Options_vec
      Options = TmbData$Options_list$Options
    }
    
    # Augment objects in TmbData (to deal with backwards compatibility)
    if( !("n_e" %in% names(TmbData)) ){
      TmbData[["n_e"]] = TmbData$n_c
    }
    if( !("ObsModel_ez" %in% names(TmbData)) ){
      TmbData[["ObsModel_ez"]] = rep(1,TmbData[["n_e"]]) %o% TmbData$ObsModel
    }
    if( !("c_iz" %in% names(TmbData)) ){
      TmbData[["c_iz"]] = matrix( TmbData$c_i, ncol=1 )
    }
    if( !("e_i" %in% names(TmbData)) ){
      TmbData[["e_i"]] = TmbData$c_iz[,1]
    }
    if( !("t_iz" %in% names(TmbData)) ){
      TmbData[["t_iz"]] = matrix( TmbData$t_i, ncol=1 )
    }
    
    # Save package version info
    capture.output( packageDescription("VAST"), file=paste0(RunDir,"/packageDescription.txt") )
    capture.output( packageDescription("FishStatsUtils"), file=paste0(RunDir,"/packageDescription.txt"), append=TRUE )
    
    # Parameters ここでエラーが出ている．make_parametersがダメ
    # TmbData=TmbData
    # if( length(Parameters)==1 && Parameters=="generate" ) Parameters = make_parameters( Version=Version, DataList=TmbData, RhoConfig=RhoConfig )
    if( length(Parameters)==1 && Parameters=="generate" ){
      Parameters = make_parameters_yk( Version=Version, DataList=TmbData, RhoConfig=RhoConfig )
      # Parameters$Beta_rho1_f = 1
      # Parameters$Beta_rho2_f = 1
      # Parameters$Epsilon_rho1_f = 1
      # Parameters$Epsilon_rho2_f = 1
    } 
    # if( length(Parameters)==1 && Parameters=="generate" ) Parameters = Return
    
    # Which parameters are turned off
    # ここでエラー
    if( length(Map)==1 && Map=="generate" ){
      Map = make_map( DataList=TmbData, TmbParams=Parameters, RhoConfig=RhoConfig, Npool=Npool )
      # Map$L_beta1_z = NULL
      # Map$L_beta2_z = NULL
      # Map$Beta_mean1_c = NULL
      # Map$Beta_mean2_c = NULL
    }else{
      warning( "Please carefully check starting values for all parameters to ensure that mapping off parameters will work as expected.")
    }
    
    # Which are random
    # Using redundant parameter names from all past versions, and then eliminating redundancies by checking against contents of Parameters
    if( length(Random)==1 && Random=="generate" ){
      Random = c("Epsiloninput1_sct", "Omegainput1_sc", "Epsiloninput1_sft", "Omegainput1_sf", "eta1_vf", "Xiinput1_scp", "Phiinput1_sk", "Epsiloninput1_sff",
                 "Epsiloninput2_sct", "Omegainput2_sc", "Epsiloninput2_sft", "Omegainput2_sf", "eta2_vf", "Xiinput2_scp", "Phiinput2_sk", "Epsiloninput2_sff",
                 "delta_i")
      if( RhoConfig[["Beta1"]]%in%c(1,2,4) ) Random = c(Random, "beta1_ct", "beta1_ft")
      if( RhoConfig[["Beta2"]]%in%c(1,2,4) ) Random = c(Random, "beta2_ct", "beta2_ft")
      if( Use_REML==TRUE ){
        Random = union(Random, c("beta1_ct","beta1_ft","gamma1_j","gamma1_tp","gamma1_ctp","lambda1_k","gamma1_cp",
                                 "beta2_ct","beta2_ft","gamma2_j","gamma2_tp","gamma2_ctp","lambda2_k","gamma2_cp"))
      }
      # Avoid problems with mapping
      Random = Random[which(Random %in% names(Parameters))]
      if( length(Random)==0) Random = NULL
    }
    
    # Bundle for debugging
    #Save = list("Map"=Map, "Data"=TmbData, "Parameters"=Parameters, "Random"=Random)
    #return(Save)
    #on.exit( return(Save) )
    #save(Save, file=paste0(RunDir,"/Save.RData"))
    
    #
    if( build_model==FALSE ){
      Return = list("Map"=Map, "Data"=TmbData, "Parameters"=Parameters, "Random"=Random)
      return( Return )
    }
    
    # Compile TMB software
    #dyn.unload( paste0(RunDir,"/",dynlib(TMB:::getUserDLL())) ) # random=Random,
    file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(CompileDir,"/",Version,".cpp"), overwrite=FALSE)
    origwd = getwd()
    on.exit(setwd(origwd),add=TRUE)
    setwd( CompileDir )
    # SEE https://github.com/kaskr/adcomp/issues/321 for flags argument
    TMB::compile( paste0(Version,".cpp"), flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
    
    # Build object
    dyn.load( paste0(CompileDir,"/",TMB::dynlib(Version)) ) # random=Random,
    Obj <- TMB::MakeADFun(data=TmbData, parameters=Parameters, hessian=FALSE, map=Map, random=Random, inner.method="newton", DLL=Version)  #
    Obj$control <- list(parscale=1, REPORT=1, reltol=1e-12, maxit=100)
    
    # Add normalization in
    if( FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v4_1_0") ){
      if( Options['normalize_GMRF_in_CPP']==FALSE ){
        message("Normalizing GMRF in R using `TMB::normalize` feature")
        Obj = TMB::normalize(Obj, flag="include_data", value=FALSE)
      }
    }
    
    # Diagnostic functions (optional)
    if( !is.null(DiagnosticDir) ){
      Obj$gr_orig = Obj$gr
      Obj$fn_orig = Obj$fn
      Obj$fn = function( vec ){
        utils::capture.output( matrix(vec,ncol=1,dimnames=list(names(Obj$par),NULL)), file=paste0(DiagnosticDir,"fn.txt") )
        utils::write.table( matrix(vec,nrow=1), row.names=FALSE, sep=",", col.names=FALSE, append=TRUE, file=paste0(DiagnosticDir,"trace.csv"))
        return( Obj$fn_orig(vec) )
      }
      Obj$gr = function( vec ){
        utils::capture.output( matrix(vec,ncol=1,dimnames=list(names(Obj$par),NULL)), file=paste0(DiagnosticDir,"gr.txt") )
        return( Obj$gr_orig(vec) )
      }
      utils::write.table( matrix(Obj$par,nrow=1), row.names=FALSE, sep=",", col.names=FALSE, file=paste0(DiagnosticDir,"trace.csv"))
    }
    
    # Local functions
    boundsifpresent_fn = function( par, map, name, lower, upper, bounds ){
      if( name %in% names(par) ){
        bounds[grep(name,names(par)),c('Lower','Upper')] = rep(1,length(grep(name,names(par)))) %o% c(lower,upper)
      }
      return( bounds )
    }
    
    # Declare upper and lower bounds for parameter search
    Bounds = matrix( NA, ncol=2, nrow=length(Obj$par), dimnames=list(names(Obj$par),c("Lower","Upper")) )
    Bounds[,'Lower'] = rep(-Inf, length(Obj$par))
    Bounds[,'Upper'] = rep( Inf, length(Obj$par))
    Bounds[grep("SigmaM",names(Obj$par)),'Upper'] = 10 # ZINB can crash if it gets > 20
    if( any(TmbData$ObsModel_ez[1,]==8) ) Bounds[grep("SigmaM",names(Obj$par)),'Upper'] = 3 # Tweedie can crash if logSigmaM gets too high
    if( !is.null(loc_x) && !is.na(Options_vec['Method']) && Options_vec['Method']==0 && Method!="Spherical_mesh" ){
      Dist = stats::dist(loc_x)
      Bounds[grep("logkappa",names(Obj$par)),'Lower'] = log( sqrt(8)/max(Dist) ) # Range = nu*sqrt(8)/kappa
      Bounds[grep("logkappa",names(Obj$par)),'Upper'] = log( sqrt(8)/min(Dist) ) # Range = nu*sqrt(8)/kappa
    }
    if( !is.na(Options_vec['Method']) && Options_vec['Method']==1 && Method!="Spherical_mesh" ){
      Bounds[grep("logkappa",names(Obj$par)),'Upper'] = log(0.9999) # Must be negative, so that Rho<1
    }
    Bounds = boundsifpresent_fn( par=Obj$par, name="ln_H_input", lower=-5, upper=5, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="gamma1", lower=-20, upper=20, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="gamma2", lower=-20, upper=20, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="lambda1", lower=-20, upper=20, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="lambda2", lower=-20, upper=20, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Beta_rho1", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Beta_rho2", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Beta_rho1_f", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Beta_rho2_f", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Epsilon_rho1", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Epsilon_rho2", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Epsilon_rho1_f", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="Epsilon_rho2_f", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="rho_c1", lower=-0.99, upper=0.99, bounds=Bounds)
    Bounds = boundsifpresent_fn( par=Obj$par, name="rho_c2", lower=-0.99, upper=0.99, bounds=Bounds)
    if( ("n_f_input"%in%names(TmbData)) && TmbData[["n_f_input"]]==0 ){
      Bounds = boundsifpresent_fn( par=Obj$par, name="L1_z", lower=c(-Inf,-0.99), upper=c(Inf,0.99), bounds=Bounds)
      Bounds = boundsifpresent_fn( par=Obj$par, name="L2_z", lower=c(-Inf,-0.99), upper=c(Inf,0.99), bounds=Bounds)
    }
    if( ("OverdispersionConfig"%in%names(TmbData)) ){
      if( TmbData[["OverdispersionConfig"]][1]==0 ) Bounds = boundsifpresent_fn( par=Obj$par, name="L1_z", lower=c(-Inf,-0.99), upper=c(Inf,0.99), bounds=Bounds)
      if( TmbData[["OverdispersionConfig"]][2]==0 ) Bounds = boundsifpresent_fn( par=Obj$par, name="L2_z", lower=c(-Inf,-0.99), upper=c(Inf,0.99), bounds=Bounds)
    }
    #for(i in 1:4){
    #  if( TmbData[["FieldConfig"]][i]==0 ){
    #    Bounds = boundsifpresent_fn( par=Obj$par, name=c("L_omega1_z","L_epsilon1_z","L_omega2_z","L_epsilon2_z")[i], lower=c(-Inf,-0.99), upper=c(Inf,0.99), bounds=Bounds)
    #  }
    #}
    
    # Change convergence tolerance
    Obj$env$inner.control$step.tol <- c(1e-8,1e-12,1e-15)[ConvergeTol] # Default : 1e-8  # Change in parameters limit inner optimization
    Obj$env$inner.control$tol10 <- c(1e-6,1e-8,1e-12)[ConvergeTol]  # Default : 1e-3     # Change in pen.like limit inner optimization
    Obj$env$inner.control$grad.tol <- c(1e-8,1e-12,1e-15)[ConvergeTol] # # Default : 1e-8  # Maximum gradient limit inner optimization
    
    # Print number of parameters
    ThorsonUtilities::list_parameters( Obj )
    
    # Return stuff
    Return = list("Obj"=Obj, "Upper"=Bounds[,'Upper'], "Lower"=Bounds[,'Lower'], "Parameters"=Parameters, "Map"=Map, "Random"=Random)
    class(Return) = "make_model"
    return( Return )
  }

#' Print model object generated by \code{\link{VAST}}
#'
#' @title Print model object
#' @param x Output from \code{\link{make_model}}
#' @param ... Not used
#' @return NULL
#' @method print make_model
#' @export
print.make_model <- function(x, ...)
{
  cat("make_model(.) result\n")
  cat("Starting values...\n")
  print( x$Obj$par )
  
  return(invisible(NULL))
}