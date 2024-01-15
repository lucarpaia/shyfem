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
  	  !! \item initialization: advertise and realize
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
        !!
	!! \subsection{Advertise}
        !!
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
	  !! the initial conditions are set and data structures are initialized:  
	  call shyfem_initialize
	  call ESMF_LogWrite("Initialized OCN", ESMF_LOGMSG_INFO, rc=rc)
      
	  ! query for importState and exportState
	  call NUOPC_ModelGet(model, importState=importState,
     +	    exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! We have added a macro to rapidly decouple the ocean component from the rest
          !! of the earth system.  Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
	  !! will result in an ocean component that does not advertise any importable
	  !! Fields. Use should you this only if you want to drive the model independently.
#define WITHIMPORTFIELDS
#ifdef WITHIMPORTFIELDS
	  !! The following are the standard variable that are imported in the 
	  !! ocean component. 
	  !! We have the atmospheric pressure at sea-level $p_a$ that
	  !! act as a forcing term on the momentum equation.
	  call NUOPC_Advertise(importState,
     +	    StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

          !! The next is the net radiation flux $R$.
	  call NUOPC_Advertise(importState,
     +	    StandardName="surface_net_downward_shortwave_flux", 
     +      name="rsns", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! We have the ocean-atmosphere flux $F_{oa}(U^{n}_o,U^{n}_a)$ 
          !! that describes the flux across the surface for momentum (wind stress) 
	  !! and temperature (heat flux). These fluxes are typically computed by the
	  !! the atmospheric component with the aid of bulk formulae.
          !! For now we have added only momentum flux that has two components.
          call NUOPC_Advertise(importState,
     +      StandardName="surface_downward_eastward_stress",
     +      name="smes", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

          call NUOPC_Advertise(importState,
     +      StandardName="surface_downward_northward_stress",
     +      name="smns", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out
#endif

	  !! As output the ocean state at the surface layer must be 
	  !! exorted. For now we have exported only the sea surface temperature.
	  call NUOPC_Advertise(exportState,
     +	    StandardName="sea_surface_temperature", name="sst", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	end subroutine

	!-----------------------------------------------------------------------------

        !!
        !! \subsection{Realize}
        !!
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
	  gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/10, 10/),
     +	minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/),
     +	maxCornerCoord=(/100000._ESMF_KIND_R8, 100000._ESMF_KIND_R8/),
     +	    coordSys=ESMF_COORDSYS_CART, 
     +      staggerLocList=(/ESMF_STAGGERLOC_CENTER/),
     +	    rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  gridOut = gridIn ! for now out same as in


          !! An ESMF Field represents a physical field, such as temperature. 
          !! The motivation for including Fields in ESMF is that bundles of Fields 
          !! are the entities that are normally exchanged when coupling Components.
          !! The ESMF Field class contains distributed and discretized field data, 
          !! a reference to its associated grid, and metadata. The Field class stores 
          !! the grid staggering for that physical field. This is the relationship 
          !! of how the data array of a field maps onto a grid (e.g. one item per cell 
          !! located at the cell center, one item per cell located at the NW corner, 
          !! one item per cell vertex, etc.). This means that different Fields which 
          !! are on the same underlying ESMF Grid but have different staggerings can 
          !! share the same Grid object without needing to replicate it multiple times.
          !! Fields can be added to States for use in inter-Component data 
          !! communications. Fields can also be added to FieldBundles, which are 
          !! groups of Fields on the same underlying Grid. One motivation for packing 
          !! Fields into FieldBundles is convenience; another is the ability to perform 
          !! optimized collective data transfers.
          !! Field communication capabilities include: data redistribution, regridding, 
          !! scatter, gather, sparse-matrix multiplication, and halo update. These are 
          !! discussed in more detail in the documentation for the specific method calls. 
          !! ESMF does not currently support vector fields, so the components of a 
          !! vector field must be stored as separate Field objects.

	  !! With the following commands we create a field for the pressure. Fields
	  !! in ESMF are of type \textsf{ESMF\_field}.
	  !! The subroutine \textsf{ESMF\_FieldCreate} simply associates the data 
	  !! with the Grid. The keywords are very important. The keyword \textsf{staggerloc}
	  !! specifies how the data is attached to the grid. \textsf{typekind} tells about 
          !! the data type, in this case double precision. The name that have been advertised must be 
	  !! also appended. 
#ifdef WITHIMPORTFIELDS

          !! An \textsf{ESMF\_Field} is created by passing the field name 
          !! (should be the same as advertised), the grid, and the data type of the 
          !! field to \textsf{ESMF\_FieldCreate}. 
	  field = ESMF_FieldCreate(name="pmsl", grid=gridIn,
     +	    typekind=ESMF_TYPEKIND_R8, staggerloc=ESMF_STAGGERLOC_CENTER, 
     +      rc=rc)
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
     +	    return  ! bail out

	  !! The others fields follow identically:
	  field = ESMF_FieldCreate(name="rsns", grid=gridIn,
     +	    typekind=ESMF_TYPEKIND_R8, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
	  call NUOPC_Realize(importState, field=field, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	  return  ! bail out

          field = ESMF_FieldCreate(name="smes", grid=gridIn,
     +      typekind=ESMF_TYPEKIND_R8, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out
          call NUOPC_Realize(importState, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +    return  ! bail out

          field = ESMF_FieldCreate(name="smns", grid=gridIn,
     +      typekind=ESMF_TYPEKIND_R8, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out
          call NUOPC_Realize(importState, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +    return  ! bail out
#endif

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

        !!
        !! \subsection{Advance SHYFEM}
        !!
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
          double precision            :: timeStepSec

!! Work in progress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	  type(ESMF_Field)            :: field
	  double precision, pointer   :: pmslPtr(:,:)
          double precision, pointer   :: smesPtr(:,:)
          double precision, pointer   :: smnsPtr(:,:)
	  integer                     :: totalLBnd(2)
	  integer                     :: totalUBnd(2)
	  integer                     :: total_count(2)
	  integer                     :: i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!! WORK IN PROGRESS

	  call ESMF_StatePrint(importState, rc=rc)

	  call ESMF_StateGet(importState, field=field, 
     +      itemName="pmsl", rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out
          
	  call ESMF_FieldGet(field, localDe=0, farrayPtr=pmslPtr,
     +	    totalLBound=totalLBnd, totalUBound=totalUBnd,
     +	    totalCount=total_count, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

          call ESMF_StateGet(importState, field=field, 
     +      itemName="smes", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out
          
          call ESMF_FieldGet(field, localDe=0, farrayPtr=smesPtr,
     +      totalLBound=totalLBnd, totalUBound=totalUBnd,
     +      totalCount=total_count, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

          call ESMF_StateGet(importState, field=field,
     +      itemName="smns", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

          call ESMF_FieldGet(field, localDe=0, farrayPtr=smnsPtr,
     +      totalLBound=totalLBnd, totalUBound=totalUBnd,
     +      totalCount=total_count, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

	  print *, "OCN: pmsl -> ", pmslPtr
          print *, "OCN: smes -> ", smesPtr
          print *, "OCN: smns -> ", smnsPtr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  !! We access to the clock object in different ways: here we
	  !! output the current time to the log file:
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

          !! Later we retrieve two important members of the clock object,
	  !! the \textsf{currTime} object and the \textsf{timeStep} object
	  call ESMF_ClockGet(clock, currTime=currTime,
     +      timeStep=timeStep, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! SHYFEM advance of one time step $\Delta t_{ao}$, we need
	  !! to access this information. From the \textsf{timeStep} 
	  !! object we query for the time step. The argument  
	  !! \textsf{s\_r8} tells to return it in second with double 
	  !! precision, which is the correct format for time variables 
          !! in SHYFEM.
	  call ESMF_TimeIntervalGet(timeStep, s_r8=timeStepSec, rc=rc)

          !! This is the call to the SHYFEM subroutine that timesteps
          !! the ocean variables for one ocean-atmosphere timestep.
	  call shyfem_run(timeStepSec)

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

        !!
        !! \subsection{Finalize SHYFEM}
        !!
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
