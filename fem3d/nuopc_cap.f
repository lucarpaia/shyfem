!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2023, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

	!! \section{The ocean component: SHYFEM}
	!!
	!! We are cascading into the objects hierarchy. We miss the last step: we need
	!! to specify the methods for the children object. Here we do this for the 
	!! ocean model object. The ocean component is coded into a Fortran module.
        module ocean_shyfem

	!! We call the modules of the libraries. First the SHYFEM library:
	use mod_shyfem
	!! and of course ESMF and NUOPC:
	use ESMF
        use NUOPC
	use NUOPC_Model, 
     +	  modelSS => SetServices

	implicit none

	private

	public SetServices

	!-----------------------------------------------------------------------------
	contains
	!-----------------------------------------------------------------------------

        !! We have already seen how to create a component with the aid of the NUOPC
        !! layer.
	subroutine SetServices(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc
        
	  rc = ESMF_SUCCESS

	  !! We do not comment again the derive/specialize commands, being
	  !! exactly equal to what we have seen for the earth component.

	  ! derive from NUOPC_Model
	  call NUOPC_CompDerive(model, modelSS, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! Instead we note that different methods are registered each one with a 
	  !! specific \textsf{specLabel} that denote a specific operation. The labels
	  !! can be regrouped into three phases: 
	  !! \begin{itemize}
  	  !! \item initialization: advertise, realize and set the clock 
	  !! \item run: advance
	  !! \item finalization: finalize
	  !! \end{itemize}

	  ! specialize model
	  call NUOPC_CompSpecialize(model, specLabel=label_Advertise,
     +	    specRoutine=Advertise, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out
	  call NUOPC_CompSpecialize(model,specLabel=label_RealizeProvided,
     +      specRoutine=Realize, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_CompSpecialize(model, specLabel=label_SetClock,
     +	    specRoutine=SetClock, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_CompSpecialize(model, specLabel=label_Advance,
     +	    specRoutine=Advance, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_CompSpecialize(model, specLabel=label_Finalize,
     +	    specRoutine=Finalize, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__)) 
     +	    return  ! bail out

	end subroutine

	!-----------------------------------------------------------------------------

	!! As done for the earth system component, the following subroutines register
	!! the members and methods for the the ocean model object that enters as first
	!! argument.
	!! The purpose of \textsf{Advertise()} is for your model to advertise its 
	!! import and export fields. This means that your model announces which model 
	!! variables it is capable of exporting (e.g., an ocean might export water 
	!! state at sea level) and which model variables it requires (e.g., an ocean 
	!! might require mass, momentum and heat fluxes as a boundary condition). The 
	!! reason there is an explicit advertise phase is because NUOPC dynamically 
	!! matches fields among all the models participating in a coupled simulation 
	!! during runtime. So, we need to collect the list of possible input and output 
	!! fields from all the models during their initialization. 
	!!  
	!! Advertising a Field does NOT allocate memory. Note that NUOPC does not 
	!! allocate memory for fields during the advertise phase or when 
	!! \textsf{NUOPC\_Advertise} is called. Instead, this is simply a way for 
	!! models to communicate the standard names of fields. During a later phase, 
	!! only those fields that are connected (e.g., a field exported from one model 
	!! that is imported by another) need to have memory allocated. Also, since ESMF
	!! will accept pointers to pre-allocated memory, it is usually not necessary to 
	!! change how memory is allocated for your model's variables. 
	subroutine Advertise(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

	  ! local variables
	  type(ESMF_State)        :: importState, exportState

	  rc = ESMF_SUCCESS

	  !! Since we are in the initialization phase, we use this subroutine 
	  !! to call also the SHYFEM initialization. Here the grid and data are read,
	  !! the initial conditions are set and data structure are initialized:  
	  call shyfem_initialize

	  call ESMF_LogWrite("Initialized OCN", ESMF_LOGMSG_INFO, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
      
	  ! query for importState and exportState
	  call NUOPC_ModelGet(model, importState=importState,
     +	    exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  ! Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
	  ! will result in a model component that does not advertise any importable
	  ! Fields. Use this if you want to drive the model independently.
#define WITHIMPORTFIELDS
#ifdef WITHIMPORTFIELDS
	  ! importable field: air_pressure_at_sea_level
	  call NUOPC_Advertise(importState,
     +	    StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  ! importable field: surface_net_downward_shortwave_flux
	  call NUOPC_Advertise(importState,
     +	    StandardName="surface_net_downward_shortwave_flux", 
     +      name="rsns", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
#endif

	  ! exportable field: sea_surface_temperature
	  call NUOPC_Advertise(exportState,
     +	    StandardName="sea_surface_temperature", name="sst", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	end subroutine

	!-----------------------------------------------------------------------------

	!! The following code fragment shows the \textsf{Realize} subroutine. During 
	!! this phase, fields that were previously advertised should now be realized. 
	!! Realizing a field means that an \textsf{ESMF\_Field} object is created and 
	!! it is added to the appropriate \textsf{ESMF\_State}, either import or export.
	subroutine Realize(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

	  ! local variables
	  type(ESMF_State)        :: importState, exportState
	  type(ESMF_TimeInterval) :: stabilityTimeStep
	  type(ESMF_Field)        :: field
	  type(ESMF_Grid)         :: gridIn
	  type(ESMF_Grid)         :: gridOut

	  rc = ESMF_SUCCESS

	  ! query for importState and exportState
	  call NUOPC_ModelGet(model, importState=importState,
     +      exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
        
	  !! In order to create an \textsf{ESMF\_Field}, you'll first need to create one 
	  !! of the ESMF geometric types, \textsf{ESMF\_Grid}, \textsf{ESMF\_Mesh}, or 
	  !! \textsf{ESMF\_LocStream}. For 2D and 3D logically rectangular grids (such as
	  !! a lat-lon grid), the typical choice is \textsf{ESMF\_Grid}. For unstructured
	  !! grids, use an \textsf{ESMF\_Mesh}. 
	  ! create a Grid object for Fields
	  gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/100, 10/),
     +	    minCornerCoord=(/10._ESMF_KIND_R8, 20._ESMF_KIND_R8/),
     +	    maxCornerCoord=(/100._ESMF_KIND_R8, 200._ESMF_KIND_R8/),
     +	    coordSys=ESMF_COORDSYS_CART, 
     +      staggerLocList=(/ESMF_STAGGERLOC_CENTER/),
     +	    rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  gridOut = gridIn ! for now out same as in

#ifdef WITHIMPORTFIELDS
	  ! importable field: air_pressure_at_sea_level
	  field = ESMF_FieldCreate(name="pmsl", grid=gridIn,
     +	    typekind=ESMF_TYPEKIND_R8, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_Realize(importState, field=field, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! An \textsf{ESMF\_Field} is created by passing the field name 
	  !! (should be the same as advertised), the grid, and the data type of the 
	  !! field to \textsf{ESMF\_FieldCreate}. 
	  ! importable field: surface_net_downward_shortwave_flux
	  field = ESMF_FieldCreate(name="rsns", grid=gridIn,
     +	    typekind=ESMF_TYPEKIND_R8, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  !! Fields are put into import or export States by calling 
	  !! \textsf{NUOPC\_Realize}.
	  call NUOPC_Realize(importState, field=field, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	  return  ! bail out
#endif

	  ! exportable field: sea_surface_temperature
	  field = ESMF_FieldCreate(name="sst", grid=gridOut,
     +	    typekind=ESMF_TYPEKIND_R8, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_Realize(exportState, field=field, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +     file=__FILE__))
     +     return  ! bail out

	end subroutine

	!-----------------------------------------------------------------------------

	subroutine SetClock(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

	  ! local variables
	  type(ESMF_Clock)              :: clock
	  type(ESMF_TimeInterval)       :: stabilityTimeStep

	  rc = ESMF_SUCCESS

	  ! query for clock
	  call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  ! initialize internal clock
	  ! here: parent Clock and stability timeStep determine actual model timeStep
	  !TODO: stabilityTimeStep should be read in from configuation
	  !TODO: or computed from internal Grid information
	  call ESMF_TimeIntervalSet(stabilityTimeStep, m=5, rc=rc) ! 5 minute steps
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_CompSetClock(model, clock, stabilityTimeStep, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	end subroutine

	!-----------------------------------------------------------------------------

	!! The SHYFEM ocean model advances the shallow water multilayer equations for
	!! stratified flows of one time step:
	!! \[
	!! U^{n+1}_o = U^{n+1}_o + \Delta t_{ao}\left( L_o(U^{n+1}_o)
	!!           + F_{oa}(U^{n}_o,U^{n}_a) \right)
	!! \]
	!! where $U_o(x,t)$ is the ocean state, $\Delta t_{ao}$ is the ocean-atmosphere
	!! timestep. In multilayer models $U_o$ pile-up the ocean state of each layer
	!! as:
	!! \[
	!! U_o = \{\, \zeta,\,U_{o,1},\,...\,U_{o,\alpha},\,...U_{o,N} \,\}
	!! \]	
	!! with $\zeta(x,t)$  the free-surface, $U_{o,\alpha}$ the momentum, 
	!! temperature and salinity of layer $\alpha$. $N$ the number of layers.
	!! The operator $L_o$ is a finite element discretization of the shallow water 
	!! multilayer equations. $F_{oa}$ are the atmosphere-ocean fluxes.
	subroutine Advance(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

  	  !! In the ESMF vocabulary the import state are atmosphere-ocean
	  !! fluxes $F_{oa}$, the export state is the oceanic state at the 
	  !! first layer $U_{o,1}$.
	  ! local variables
	  type(ESMF_Clock)            :: clock
    	  type(ESMF_State)            :: importState, exportState
	  type(ESMF_Time)             :: currTime
	  type(ESMF_TimeInterval)     :: timeStep
	  character(len=160)          :: msgString
	  double precision, parameter :: minusOne = -1.

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionEnter("OCN:Advance")
#endif

	  rc = ESMF_SUCCESS

	  ! query for clock, importState and exportState
	  call NUOPC_ModelGet(model, modelClock=clock,
     +      importState=importState, exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! This is the call to the SHYFEM subroutine that timesteps
	  !! the ocean variables for one ocean atmosphere timestep.
	  call shyfem_run(minusOne)

	  ! Because of the way that the internal Clock was set in SetClock(),
	  ! its timeStep is likely smaller than the parent timeStep. As a consequence
	  ! the time interval covered by a single parent timeStep will result in
	  ! multiple calls to the Advance() routine. Every time the currTime
	  ! will come in by one internal timeStep advanced. This goes until the
	  ! stopTime of the internal Clock has been reached.

     	  call ESMF_ClockPrint(clock, options="currTime",
     +	    preString="------>Advancing OCN from: ",unit=msgString,rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  call ESMF_ClockGet(clock, currTime=currTime,
     +      timeStep=timeStep, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  call ESMF_TimePrint(currTime + timeStep,
     +	    preString="---------------------> to: ",unit=msgString,rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
 	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionExit("OCN:Advance")
#endif
	end subroutine

	!-----------------------------------------------------------------------------

  	!! We register a method for the finalization of the ocean code. Here files are 
  	!! closed and the processes is killed with the memory cleaned through the 
	!! deallocations of the data structures. 
	subroutine Finalize(model, rc)
	  type(ESMF_GridComp)  :: model
	  integer, intent(out) :: rc

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionEnter("OCN:Finalize")
#endif

	  rc = ESMF_SUCCESS

	  ! HERE THE MODEL IS FINALIZED: 
	  call shyfem_finalize

	  call ESMF_LogWrite("Finalized OCN", ESMF_LOGMSG_INFO, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +	    return  ! bail out

#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionExit("OCN:Finalize")
#endif
	end subroutine

	!-----------------------------------------------------------------------------

	end module
