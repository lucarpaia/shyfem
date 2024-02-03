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

	!! We declare private members and methods. \textsf{SHYFEM\_FieldMetadata}
	!! is an object that was not present in ESMF. We have created it in the cap
	!! layer and for this reason it has the SHYFEM prefix. This object contains
	!! important information to create the ESMF fields and will help us
	!! in factorizing field creation within do loops, instead of doing evreything
	!! manually. The info concerns the field names and mesh staggering
	!! and are initialized here with default values.
	private

	type SHYFEM_Metadata
	  character(4)                :: fieldName = " "
          character(50)               :: longName = " "
	  type(ESMF_MeshLoc)          :: meshloc = ESMF_MESHLOC_NODE
	end type

	integer :: SHYFEM_numOfImportFields
        integer :: SHYFEM_numOfExportFields
	type(SHYFEM_Metadata), allocatable, 
     +     save :: SHYFEM_FieldMetadata(:)

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
	  integer                 :: var, numImp

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

	  !! The fields are organized with a precise order in a \textsf{SHYFEM\_FieldMetadata}
	  !! array. This ordering is beneficial to factorize the code and perform the 
	  !! many operations of the fields with do loops. You have to keep in mind that, 
	  !! first we have import fields, and only after we have export fields. 
	  allocate( SHYFEM_FieldMetadata(8) )

          !! We have added two macros to rapidly decouple the ocean component from the rest
          !! of the earth system. Disabling both the following macros,
          !! will result in an ocean component that does not advertise any importable
          !! Fields. Use should you this only if you want to drive the model independently.
          !! You could also couples the ocean component through momentum or heat flux. 
          !! This may be useful for debugging, for exmaple. In general, both the macros must
          !! be defined.
#define WITHIMPORTFIELDS_MOMENTUMFLUX
#define WITHIMPORTFIELDS_HEATFLUX

          !! The metadata of each field is filled with a simple constructor statement.
          !! In the following we describe each field rapidly. 
          numImp = 0

#ifdef WITHIMPORTFIELDS_MOMENTUMFLUX
          numImp = 3

          !! \begin{itemize}
          !! \item atmospheric pressure at sea-level $p_a$ that
          !! act as a forcing term on the momentum equation.
          SHYFEM_FieldMetadata(1) = SHYFEM_Metadata(fieldName="pmsl",
     +      longName="air_pressure_at_sea_level",
     +      meshloc=ESMF_MESHLOC_NODE)

          !! \item ocean-atmosphere flux $F_{oa}(U^{n}_o,U^{n}_a)$ 
          !! that describes the flux across the surface for momentum (wind stress) 
          !! and temperature (heat flux). These fluxes are typically computed by the
          !! the atmospheric component with the aid of bulk formulae.
          !! For now we have added only momentum flux that has two components.
          SHYFEM_FieldMetadata(2) = SHYFEM_Metadata(fieldName="smes",
     +      longName="surface_downward_eastward_stress",
     +      meshloc=ESMF_MESHLOC_NODE)

          SHYFEM_FieldMetadata(3) = SHYFEM_Metadata(fieldName="smns",
     +      longName="surface_downward_northward_stress",
     +      meshloc=ESMF_MESHLOC_NODE)

#endif
#ifdef WITHIMPORTFIELDS_HEATFLUX
          numImp = numImp + 1
          !! \item net radiation flux $R$.
          SHYFEM_FieldMetadata(numImp) = SHYFEM_Metadata(
     +      fieldName="rsns",
     +      longName="surface_net_downward_shortwave_flux",
     +      meshloc=ESMF_MESHLOC_NODE)

	  numImp = numImp + 1
          !! \item sensible heat flux $R_{long}$.
          SHYFEM_FieldMetadata(numImp) = SHYFEM_Metadata(
     +      fieldName="rlns",
     +      longName="surface_net_downward_longwave_flux",
     +      meshloc=ESMF_MESHLOC_NODE)

          !! \item sensible heat flux $Q_{sens}$.
!          SHYFEM_FieldMetadata(numImp) = SHYFEM_Metadata(
!     +      fieldName="stsh",
!     +      longName="surface_downward_sensible_heat_flux_in_air",
!     +      meshloc=ESMF_MESHLOC_NODE)

          !! \item sensible heat flux $Q_{lat}$.
!          SHYFEM_FieldMetadata(numImp) = SHYFEM_Metadata(
!     +      fieldName="stlh",
!     +      longName="surface_downward_heat_flux_in_air",
!     +      meshloc=ESMF_MESHLOC_NODE)
#endif
          SHYFEM_numOfImportFields = numImp
          SHYFEM_numOfExportFields = 1

	  !! \item sea surface salinity.
	  SHYFEM_FieldMetadata(numImp+1) = SHYFEM_Metadata(
     +      fieldName="sst",
     +      longName="sea_surface_temperature",
     +      meshloc=ESMF_MESHLOC_NODE)

          !! \end{itemize}
	  !! The with a do loop The following we advertise the import state in the 
	  !! ocean component.
	  do var=1,SHYFEM_numOfImportFields

	    call NUOPC_Advertise(importState,
     +	      StandardName=SHYFEM_FieldMetadata(var)%longName, 
     +        name=SHYFEM_FieldMetadata(var)%fieldName, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out

	  enddo

	  !! As output the ocean state at the surface layer must be 
	  !! exorted. For now we have exported only the sea surface temperature.
          do var=SHYFEM_numOfImportFields+1,
     +         SHYFEM_numOfImportFields+SHYFEM_numOfExportFields

	    call NUOPC_Advertise(exportState,
     +        StandardName=SHYFEM_FieldMetadata(var)%longName,
     +        name=SHYFEM_FieldMetadata(var)%fieldName, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out

	  enddo

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
!	  type(ESMF_Grid)         :: gridIn
!	  type(ESMF_Grid)         :: gridOut
          type(ESMF_Mesh)         :: meshIn
	  integer                 :: var

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
	  !! SHYFEM has already read the mesh during the advertise phase. With the command
	  !! \textsf{SHYFEM\_MeshGet} the mesh is imported into a \textsf{ESMF\_Mesh} 
	  !! object and later we also write into vtk file to check that evreything is ok.
	  !! Vtk files can be opened with Paraview.
!	  gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/10, 10/),
!     +	minCornerCoord=(/-5._ESMF_KIND_R8, -5._ESMF_KIND_R8/),
!     +	maxCornerCoord=(/5._ESMF_KIND_R8, 5._ESMF_KIND_R8/),
!     +	    coordSys=ESMF_COORDSYS_CART, 
!     +      staggerLocList=(/ESMF_STAGGERLOC_CENTER/),
!     +	    rc=rc)
!	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
!     +	    line=__LINE__,
!     +	    file=__FILE__))
!     +	    return  ! bail out
!	  gridOut = gridIn ! for now out same as in

          call SHYFEM_MeshGet(meshIn, rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return
          call ESMF_LogWrite("  Get Mesh OCN", ESMF_LOGMSG_INFO, rc=rc)

          call ESMF_MeshWrite(meshIn, "shyfem_mesh", rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

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
	  !!
	  !! With the following commands we create a field for the pressure. Fields
	  !! in ESMF are of type \textsf{ESMF\_field}.
	  !! The subroutine \textsf{ESMF\_FieldCreate} simply associates the data 
	  !! with the Grid. The keywords are very important. The keyword \textsf{staggerloc}
	  !! specifies how the data is attached to the grid. \textsf{typekind} tells about 
          !! the data type, in this case double precision. The name that have been advertised must be 
	  !! also appended.
	  !! 
          !! An \textsf{ESMF\_Field} is created by passing the field name 
          !! (should be the same as advertised), the grid, the staggering, and the data type of the 
          !! field to \textsf{ESMF\_FieldCreate}.
	  do var=1,SHYFEM_numOfImportFields

	    field = ESMF_FieldCreate(meshIn,
     +	      ESMF_TYPEKIND_R8,
     +        meshloc=SHYFEM_FieldMetadata(var)%meshloc, 
     +        name=SHYFEM_FieldMetadata(var)%fieldName,
     +        rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out

            !! Fields are put into import or export States by calling 
            !! \textsf{NUOPC\_Realize}.
	    call NUOPC_Realize(importState, field=field, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out

	  enddo

          do var=SHYFEM_numOfImportFields+1,
     +         SHYFEM_numOfImportFields+SHYFEM_numOfExportFields

	    field = ESMF_FieldCreate(meshIn,
     +        ESMF_TYPEKIND_R8,
     +        meshloc=SHYFEM_FieldMetadata(var)%meshloc,
     +        name=SHYFEM_FieldMetadata(var)%fieldName,
     +        rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out

	    call NUOPC_Realize(exportState, field=field, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +       file=__FILE__))
     +       return  ! bail out

	  enddo

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
	!! In this subrotine we do two things. We take the set the fluxes $F_{oa}$
	!! inside SHYFEM, through a quite complex chain of commands. Then we advance 
	!! SHYFEM.
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

	  ! local variables 
	  type(ESMF_Field)            :: field
	  double precision, pointer   :: fieldPtr(:)
	  integer                     :: var

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
	  call ESMF_TraceRegionEnter("OCN:Advance")
#endif

	  rc = ESMF_SUCCESS

	  !! The various objects are nested. We access the fields with three
	  !! consecutive \textsf{Get} command. We start accessing the most
	  !! exterior component, namely the model. We query for three
	  !! objects: the \textsf{clock}, the \textsf{importState} 
	  !! and \textsf{exportState}.
	  call NUOPC_ModelGet(model, modelClock=clock,
     +      importState=importState, exportState=exportState, rc=rc)
	  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	    line=__LINE__,
     +	    file=__FILE__))
     +	    return  ! bail out

	  !! With a little bit of paranoia we print the import state object
	  !! (this needs to be removed or added only in a debug mode):
	  call ESMF_StatePrint(importState, rc=rc)

	  !! From the component, we retrieve the state and from the state
	  !! we retreive the field. The field data is typically a large array 
	  !! stored in memory, we use a pointer to it. Field can be written to a file.
	  !! Finally, after all these \textsf{Get} commands, a \textsf{Set} commands
	  !! copy the field pointer to the SHYFEM flux variables.
	  do var=1,SHYFEM_numOfImportFields

	    call ESMF_StateGet(importState, field=field, 
     +        itemName=SHYFEM_FieldMetadata(var)%fieldName, rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out
          
	    call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr,
     +	      rc=rc)
	    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +	      line=__LINE__,
     +	      file=__FILE__))
     +	      return  ! bail out

	    call SHYFEM_FieldWrite(field, 
     +      trim(SHYFEM_FieldMetadata(var)%fieldName)//trim(".vtk"), rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +        line=__LINE__,
     +        file=__FILE__))
     +        return  ! bail out

	    print *, fieldPtr

	    call SHYFEM_FieldSet(fieldPtr,
     +        SHYFEM_FieldMetadata(var)%fieldName, rc)
            if(ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,
     +        line=__LINE__,
     +        file=__FILE__))
     +        return  ! bail out

	  enddo

	  !! The first task is completed. We move to the second one:
	  !! advancing the model of one timestep.
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
     +      file=__FILE__))
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

	  deallocate(SHYFEM_FieldMetadata)

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

        !!
        !! \subsection{Get SHYFEM mesh}
        !!
	!! This section describes how to get the mesh from SHYFEM and set into a ESMF 
	!! Mesh class. It starts with an explanation of creating a Mesh and then goes 
	!! through other Mesh methods. To create a Mesh we need to set some properties 
	!! of the Mesh as a whole, some properties of each node in the mesh and then 
	!! some properties of each element which connects the nodes.
	!! For the Mesh as a whole we set its parametric dimension (\textsf{parametricDim}) 
	!! and spatial dimension (\textsf{spatialDim}). A Mesh's parametric dimension can 
	!! be thought of as the dimension of the elements which make up the Mesh.
	!! A Mesh's spatial dimension, on the other hand, is the is the number of coordinate 
	!! dimensions needed to describe the location of the nodes making up the Mesh. For
	!! SHYFEM both are always two, which means that the mesh is defined in parametric
	!! two-dimensional manifold. Even is spherical coordinates, SHYFEM use the spherical
	!! transformation to transform the differential operator on a two-dimensional. 
	!! Another important properties of SHYFEM mesh is that they are made of only 
	!! triangular elements. Without loss of generality we can assume that a mesh looks
	!! like the one in figure:
	!!
	!!\begin{minipage}{\linewidth} 
	!!\begin{verbatim}
	!!
	!! 
	!!  2.0   7 ------- 8 ------- 9
	!!        |  \   6  |  \   4  |
	!!        |    \    |    \    |
	!!        |  7    \ |  5    \ |
	!!  1.0   4 ------- 5 ------- 6
	!!        |  \   8  |  \   3  |
	!!        |    \    |    \    |
	!!        |  1    \ |  2   \  |
	!!  0.0   1 ------- 2 ------- 3
	!!
	!!       0.0       1.0        2.0 
	!! 
	!!        Node Id labels at corners
	!!       Element Id labels in centers
	!!       (Everything owned by PET 0) 
	!!
	!!\end{verbatim}
	!!\end{minipage}
	!!
	!! With this is mind let's have a look to the subroutine that get SHYFEM mesh.
	subroutine SHYFEM_MeshGet(SHYFEM_mesh, rc)

          type(ESMF_Mesh), intent(out)    :: SHYFEM_mesh
          integer, intent(out) 		  :: rc

	  !! First note the type \textsf{ESMF\_Mesh} which defines an unstructured grid
	  !! object. Other global and local properties of the mesh follows: 
          integer                         :: SHYFEM_spatialDim = 2
	  integer                         :: SHYFEM_numberOfNodes
	  integer                         :: SHYFEM_numberOfElements
	  integer, allocatable            :: SHYFEM_nodeIds(:)
	  integer, allocatable            :: SHYFEM_nodeOwners(:)
	  real(ESMF_KIND_R8), allocatable :: SHYFEM_nodeCoords(:)
          integer, allocatable            :: SHYFEM_elementIds(:)
	  integer, allocatable            :: SHYFEM_elementTypes(:)
	  integer, allocatable            :: SHYFEM_elementConn(:)

	  integer :: i, ie

	  rc = ESMF_SUCCESS

	  !! We save the number of nodes and the number of elements
	  !! We got this from the well known global variables \textsf{nkn}
	  !! and \textsf{nel} of SHYFEM. These are available into the module
	  !! \textsf{grid} which is nested into \textsf{mod\_shyfem}
	  SHYFEM_numberOfNodes    = nkn
          SHYFEM_numberOfElements = nel

          allocate( SHYFEM_nodeIds(SHYFEM_numberOfNodes) )
	  allocate( SHYFEM_nodeCoords(SHYFEM_numberOfNodes
     +      *SHYFEM_spatialDim) )
          allocate( SHYFEM_nodeOwners(SHYFEM_numberOfNodes) )

	  allocate( SHYFEM_elementIds(SHYFEM_numberOfElements) )
          allocate( SHYFEM_elementTypes(SHYFEM_numberOfElements) )
          allocate( SHYFEM_elementConn(SHYFEM_numberOfElements*3) )

	  !! The structure of the per node and element information used to 
	  !! create a Mesh is influenced by the Mesh distribution strategy. 
	  !! The Mesh class is distributed by elements. This means that a node 
	  !! must be present on any PET that contains an element associated with 
	  !! that node, but not on any other PET (a node can't be on a PET
	  !! without an element "home"). Since a node may be used by two or more 
	  !! elements located on different PETs, a node may be duplicated on 
	  !! multiple PETs. When a node is duplicated in this manner, 
	  !! one and only one of the PETs that contain the node must "own" the node. 
	  !! The user sets this ownership when they define the nodes during Mesh 
	  !! creation. For each node in the Mesh we set three properties:
	  !! the global id of the node (\textsf{nodeIds}), node coordinates 
	  !! (\textsf{nodeCoords}), and which PET owns the node ({\textsf{nodeOwners}).
	  !! The node id is a unique (across all PETs) integer attached to 
	  !! the particular node. It is used to indicate which nodes are the 
	  !! same when connecting together pieces of the Mesh on different 
	  !! processors. The node coordinates indicate the location of a node 
	  !! in space and are used in the \textsf{ESMF\_FieldRegrid()} functionality 
	  !! when interpolating. The node owner indicates which PET is in charge 
	  !! of the node. This is used when creating a Field on the Mesh to 
	  !! indicate which PET should contain a Field location for the data. 
	  !! The \textsf{nodeCoords} is an array containing the physical coordinates 
	  !! of the nodes to be created on this PET. This input consists of a 1D array 
	  !! the size of the number of nodes on this PET times the Mesh's spatial 
	  !! dimension (\textsf{spatialDim}). The coordinates in this array are ordered 
	  !! so that the coordinates for a node lie in sequence in memory. For the 
	  !! example shown above coordinates are in the following sequence
	  !! $\{x_1,\,y_1,\,...,\,x_9,\,y_9\}$. 
	  SHYFEM_nodeOwners = 0
	  do i=1,SHYFEM_numberOfNodes
	    SHYFEM_nodeIds(i) = i
	    SHYFEM_nodeCoords(i*SHYFEM_spatialDim-1) = xgv(i)
	    SHYFEM_nodeCoords(i*SHYFEM_spatialDim)   = ygv(i)
	  enddo

	  !! For each element in the Mesh we set three properties: the global id 
	  !! of the element (\textsf{elementIds}), the topology type of the element 
	  !! (\textsf{elementTypes}), and which nodes are connected together to form 
	  !! the element (\textsf{elementConn}). The element id is a unique (across all PETs) 
	  !! integer attached to the particular element. The element type describes 
	  !! the topology of the element (e.g. a triangle vs. a quadrilateral). 
	  !! The range of choices for the topology of the elements in a Mesh are 
	  !! The element connectivity indicates which nodes are to be connected 
	  !! together to form the element. The number of nodes connected together for 
	  !! each element is implied by the elements topology type ({\textsf{elementTypes}). 
	  !! It is IMPORTANT to note, that 
	  !! the entries in this list are NOT the global ids of the nodes, but are indices 
	  !! into the PET local lists of node info used in the Mesh Create. In other words, 
	  !! the element connectivity isn't specified in terms of the global list of nodes, 
	  !! but instead is specified in terms of the locally described node info. One 
	  !! other important point about connectivities is that the order of the nodes in the 
	  !! connectivity list of an element is important. In general, when specifying an 
	  !! element with parametric dimension 2, the nodes should be given in 
	  !! counterclockwise order around the element. 
	  SHYFEM_elementTypes = ESMF_MESHELEMTYPE_TRI
          do ie=1,SHYFEM_numberOfElements
            SHYFEM_elementIds(ie) = ie
	    do i=1,3
              SHYFEM_elementConn((ie-1)*3+i) = nen3v(i,ie)
	    enddo
          enddo

	  !! Once we have collected the mesh properties in the correct form, 
	  !! a final call create the mesh object at once. We also print to
	  !! file in vtk format, that can be visualize with Paraview.
	  SHYFEM_mesh = ESMF_MeshCreate(parametricDim=SHYFEM_spatialDim,
     +      spatialDim=SHYFEM_spatialDim,
     +      coordSys=ESMF_COORDSYS_CART,
     +      nodeIds=SHYFEM_nodeIds, nodeCoords=SHYFEM_nodeCoords,
     +      nodeOwners=SHYFEM_nodeOwners, elementIds=SHYFEM_elementIds,
     +      elementTypes=SHYFEM_elementTypes, 
     +      elementConn=SHYFEM_elementConn,
     +      rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

          deallocate( SHYFEM_nodeIds, SHYFEM_nodeOwners, 
     +      SHYFEM_nodeCoords )
	  deallocate( SHYFEM_elementIds, SHYFEM_elementTypes, 
     +      SHYFEM_elementConn )

	end subroutine

        !-----------------------------------------------------------------------------

        !! Eventually we write the regridded fields to files.
	!! This can be helpful for debugging and checking the interpolations.
	!! We can write such a file with ESMF subroutine \textsf{SHYFEM\_FieldWrite}
	!! but this works only with the third party library PARALLELIO (PIO).
	!! Moreover the only format allowed when this manual was written was netcdf
	!! (ugrid). We have preferred to do use vtk format to visualize the data,
	!! as done with the mesh. This lead to only one type of file outputted.
	!! Vtk files can be visualize nicely with Paraview.
        subroutine SHYFEM_FieldWrite(field, filename, rc)

          type(ESMF_Field), intent(in)  :: field
	  character(len=*), intent(in)  :: filename
          integer, intent(out)          :: rc

	  ! local variables
          double precision, pointer   :: ptr(:)
          character(20) :: str
          integer :: i, ie, ii

          rc = ESMF_SUCCESS

	  !! We recovet the field with the usual call
          call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr,
     +      rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,
     +      line=__LINE__,
     +      file=__FILE__))
     +      return  ! bail out

	  !! We spend some word about the vtk format. We use
	  !! three digits to represents point coordinates, that means
	  !! that you can reach mm of resolution in you mesh but not beyond.
	  !! Paraview does not like zeros as float 0.00 and I had to manually
	  !! write zero as an integer. For the variable precision we assume
	  !! that three digits are also sufficient.
	  open(1, file = filename, status = 'unknown')
	  write(1,'(a)') "# vtk DataFile Version 3.0"
	  write(1,'(a)') "This file was generated by NUOPC"
	  write(1,'(a)') "ASCII"
	  write(1,'(a)') "DATASET UNSTRUCTURED_GRID"
	  write(1,'(a)', advance='no') "POINTS ";
	  write(1, *) nkn, " double"
	  do i=1,nkn
            write(str,'(f6.3)') xgv(i)
	    str = trim(adjustl(str))
	    do ii = len_trim(str),1,-1
	      if (str(ii:ii)/="0") exit
	    enddo
            if (str(ii:ii)==".") ii=ii-1
	    write(1,'(aX)', advance="no") str(1:ii)
            write(str,'(f6.3)') ygv(i)
            str = trim(adjustl(str))
            do ii = len_trim(str),1,-1
              if (str(ii:ii)/="0") exit
            enddo
	    if (str(ii:ii)==".") ii=ii-1
            write(1,'(aX,i1)') str(1:ii), 0
	  end do
	  write(1,'(a)', advance='no') "CELLS "
	  write(1, *) nel, nel*4
          do ie=1,nel
            write(1,"(i6X,i6X,i6X,i6)") 3,
     +	      nen3v(1,ie)-1, nen3v(2,ie)-1, nen3v(3,ie)-1
          end do
          write(1,'(a)', advance="no") "CELL_TYPES "
	  write(1, *) nel
          do ie=1,nel
            write(1,*) 5
          end do
	  write(1,'(a)', advance="no") "POINT_DATA "
	  write(1, *) nkn
	  write(1,'(a)') "SCALARS _NODE_NUM double 1"
	  write(1,'(a)') "LOOKUP_TABLE default"
          do i=1,nkn
            write(str,'(f11.4)') ptr(i)
            str = trim(adjustl(str))
            do ii = len_trim(str),1,-1
              if (str(ii:ii)/="0") exit
            enddo
            if (str(ii:ii)==".") ii=ii-1
            write(1,'(aX)', advance="no") str(1:ii)
          end do
	  close(1)

        end subroutine

        !-----------------------------------------------------------------------------

        !! Eventually we write the regridded fields to files.
        !! This can be helpful for debugging and checking the interpolations.
        !! We can write such a file with ESMF subroutine \textsf{SHYFEM\_FieldWrite}
        !! but this works only with the third party library PARALLELIO (PIO).
        !! Moreover the only format allowed when this manual was written was netcdf
        !! (ugrid). We have preferred to do use vtk format to visualize the data,
        !! as done with the mesh. This lead to only one type of file outputted.
        !! Vtk files can be visualize nicely with Paraview.
        subroutine SHYFEM_FieldSet(fieldPtr, fieldName, rc)

          double precision, pointer, intent(in) :: fieldPtr(:)
          character(len=*), intent(in)          :: fieldName
          integer, intent(out)                  :: rc

          rc = ESMF_SUCCESS

	  select case (fieldName)
	    case ("pmsl")
	      ppv = fieldPtr
	    case ("smes")
	      tauxnv = fieldPtr
            case ("smns")
	      tauynv = fieldPtr
	    case ("rsns")
	      metrad = fieldPtr
            case ("rlns")
              !metrad = fieldPtr
	    case default
              call ESMF_LogWrite("  OCN unknown field name", 
     +          ESMF_LOGMSG_INFO, rc=rc)
	      rc = 1
	  end select

	end subroutine

	end module
