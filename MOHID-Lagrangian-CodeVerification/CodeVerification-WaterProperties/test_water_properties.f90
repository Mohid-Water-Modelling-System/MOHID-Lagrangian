program test_water_properties
    use simulationPrecision_mod, only: prec
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use simulationGlobals_mod, only: Globals
    use simulationInitialize_mod, only: InitFromXml
    use kernelVerticalMotion_mod, only: kernelVerticalMotion_class
    use kernelFloater_mod, only: kernelFloater_class
    use stateVector_mod, only: stateVector_class
    use background_mod, only: background_class
    use stringifor, only: string
    use xmlparser_mod, only: XMLReader
    use FoX_dom, only: Node, getDocumentElement
    use legacy_water_reference
    implicit none
    type(kernelVerticalMotion_class) :: water
    type(kernelFloater_class) :: floater
    type(stateVector_class) :: sv, buoySV
    type(background_class) :: emptyBackground(0)
    type(string) :: text, tag, attribute, path
    type(Node), pointer :: document, constants
    real(prec) :: S(4), T(4), Z(4), rho(4), mu(4), baseRho(4), baseMu(4)
    real(prec), allocatable :: nu(:), actualRho(:)
    integer :: d, v, i, checks
    real(prec) :: referenceMu, ratio(6), expectedKinematicViscosity
    real(prec) :: velocity(4,14), usedRho(4), usedNu(4), contrast(4), diameter(4), expectedSpeed(4)
    logical :: found
    character(len=256) :: option, argument
    call get_command_argument(1,option)
    call get_command_argument(2,argument)
    ! Invalid-input modes are called in separate processes by the runner.
    if (trim(option)=='density-flag') then
        text=trim(argument); call Globals%Constants%setUseConstantWaterDensity(text); stop
    elseif (trim(option)=='viscosity-flag') then
        text=trim(argument); call Globals%Constants%setUseConstantWaterDynamicViscosity(text); stop
    elseif (trim(option)=='mean-density') then
        text=trim(argument); call Globals%Constants%setMeanDensity(text); stop
    elseif (trim(option)=='dynamic-viscosity') then
        text=trim(argument); call Globals%Constants%setMeanDynamicViscosity(text); stop
    endif
    if (trim(option)=='kinematic-xml') then
        path='./'
        call Globals%initialize(path)
        path=trim(argument)
        call InitFromXml(path)
        call get_command_argument(3,argument)
        read(argument,*) expectedKinematicViscosity
        if (abs(Globals%Constants%MeanKinematicViscosity-expectedKinematicViscosity)>1e-18_prec) &
            error stop 'Kinematic-viscosity XML alias or priority is incorrect'
        print *, 'PASS: kinematic viscosity XML compatibility'
        stop
    endif
    checks=0
    call require(.not.Globals%Constants%UseConstantWaterDensity,'density override defaults off')
    call require(.not.Globals%Constants%UseConstantWaterDynamicViscosity,'viscosity override defaults off')
    Globals%Constants%MeanDensity=1030._prec
    Globals%Constants%MeanKinematicViscosity=1.09e-6_prec
    Globals%Var%temp='temp'; Globals%Var%sal='sal'
    S=[35._prec,35._prec,30._prec,37._prec]
    T=[15._prec,15._prec,5._prec,25._prec]
    Z=[0._prec,-1000._prec,-20._prec,-200._prec]
    baseRho=ref_seaWaterDensity(S,T,Z)
    baseMu=ref_absoluteSeaWaterViscosity(S,T)
    call require(baseRho(2)>baseRho(1),'pressure changes baseline density at equal T/S')
    call require(baseMu(1)==baseMu(2),'dynamic viscosity has no depth dependence')
    allocate(sv%state(4,5),sv%varName(5))
    sv%varName(1)='x';sv%varName(2)='y';sv%varName(3)='z'
    sv%varName(4)='temp';sv%varName(5)='sal'
    sv%state=0;sv%state(:,3)=Z;sv%state(:,4)=T;sv%state(:,5)=S
    text='0.001';call Globals%Constants%setMeanDynamicViscosity(text)
    do d=0,1
        text=merge('1','0',d==1);call Globals%Constants%setUseConstantWaterDensity(text)
        do v=0,1
            text=merge('1','0',v==1);call Globals%Constants%setUseConstantWaterDynamicViscosity(text)
            rho=water%seaWaterDensity(S,T,Z)
            mu=water%absoluteSeaWaterViscosity(S,T)
            if (d==1) then
                call require(all(rho==Globals%Constants%MeanDensity),'constant density independent of T/S/Z')
            else
                call require(maxval(abs(rho-baseRho))<1e-11_prec,'variable density matches original equation')
            endif
            if (v==1) then
                call require(all(mu==Globals%Constants%MeanDynamicViscosity),'constant mu independent of T/S')
            else
                call require(maxval(abs(mu-baseMu))<1e-16_prec,'variable mu matches original equation')
            endif
            call floater%DensityKViscosity_seawater(sv,nu,actualRho)
            call require(maxval(abs(actualRho-rho))<1e-11_prec,'floater retains independently selected density')
            call require(maxval(abs(nu-mu/rho))<1e-18_prec,'kinematic viscosity is always selected mu/rho')
        enddo
    enddo
    ! A high prescribed mu is used directly, without a separate kinematic limiter.
    Globals%Constants%MeanDynamicViscosity=0.01_prec
    Globals%Constants%UseConstantWaterDensity=.false.
    Globals%Constants%UseConstantWaterDynamicViscosity=.true.
    Globals%Constants%MeanKinematicViscosity=1e-10_prec
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(maxval(abs(nu-0.01_prec/baseRho))<1e-18_prec,'MeanKinematicViscosity does not override prescribed mu/rho')
    call require(nu(2)<nu(1),'constant mu with depth-dependent density gives variable nu')
    Globals%Constants%UseConstantWaterDensity=.true.
    Globals%Constants%MeanDensity=100._prec
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(all(actualRho==100._prec),'selected density is independent of viscosity')
    call require(all(abs(nu-0.0001_prec)<1e-18_prec),'both constant properties give constant nu')
    ! Exercise low/high dynamic-viscosity fallback for both density choices.
    Globals%Constants%MeanDensity=1030._prec
    do i=0,1
        referenceMu=0.1_prec
        if(i==1)referenceMu=0.0001_prec
        Globals%Constants%MeanDynamicViscosity=referenceMu
        do d=0,1
            Globals%Constants%UseConstantWaterDensity=(d==1)
            do v=0,1
                Globals%Constants%UseConstantWaterDynamicViscosity=(v==1)
                rho=baseRho
                if(d==1)rho=1030._prec
                call floater%DensityKViscosity_seawater(sv,nu,actualRho)
                call require(maxval(abs(actualRho-rho))<1e-11_prec,'mu fallback preserves chosen density')
                call require(maxval(abs(nu-referenceMu/rho))<1e-18_prec,'guard acts on mu before division')
            enddo
        enddo
    enddo
    ! Values either side of the 0.1 and 1.9 relative-distance thresholds.
    Globals%Constants%UseConstantWaterDynamicViscosity=.false.
    ratio=[0.099_prec,0.101_prec,1._prec,1.899_prec,1.901_prec,2._prec]
    do i=1,size(ratio)
        referenceMu=baseMu(1)/ratio(i)
        Globals%Constants%MeanDynamicViscosity=referenceMu
        mu=water%absoluteSeaWaterViscosity(S,T)
        if(i==1.or.i>=5)then
            call require(abs(mu(1)-referenceMu)<1e-18_prec,'out-of-range mu uses mean')
        else
            call require(abs(mu(1)-baseMu(1))<1e-18_prec,'in-range mu is retained')
        endif
    enddo
    Globals%Constants%MeanDensity=1030._prec
    Globals%Constants%MeanKinematicViscosity=1.09e-6_prec
    Globals%Constants%MeanDynamicViscosity=0.001_prec
    ! Missing fields independently use 15 C and salinity 36; property switches still apply.
    sv%varName(4)='unused1';sv%varName(5)='unused2'
    do d=0,1
        Globals%Constants%UseConstantWaterDensity=(d==1)
        rho=ref_seaWaterDensity(0._prec*Z+36._prec,0._prec*Z+15._prec,Z)
        if(d==1)rho=Globals%Constants%MeanDensity
        do v=0,1
            Globals%Constants%UseConstantWaterDynamicViscosity=(v==1)
            Globals%Constants%MeanDynamicViscosity=0.001_prec
            mu=ref_absoluteSeaWaterViscosity(0._prec*Z+36._prec,0._prec*Z+15._prec)
            if(v==1)mu=0.001_prec
            call floater%DensityKViscosity_seawater(sv,nu,actualRho)
            call require(maxval(abs(actualRho-rho))<1e-11_prec,'missing T/S respects density switch')
            call require(maxval(abs(nu-mu/rho))<1e-18_prec,'missing T/S respects independent viscosity switch')
        enddo
    enddo
    Globals%Constants%UseConstantWaterDensity=.false.
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(actualRho(2)>actualRho(1),'default T/S retains pressure dependence with depth')
    sv%varName(4)='temp'
    rho=ref_seaWaterDensity(0._prec*Z+36._prec,T,Z)
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(maxval(abs(actualRho-rho))<1e-11_prec,'temperature retained when only salinity missing')
    sv%varName(4)='unused1';sv%varName(5)='sal'
    rho=ref_seaWaterDensity(S,0._prec*Z+15._prec,Z)
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(maxval(abs(actualRho-rho))<1e-11_prec,'salinity retained when only temperature missing')
    Globals%Constants%MeanDynamicViscosity=0.001_prec
    ! Computed NaN/infinite-temperature viscosity falls back to the reference mu.
    sv%varName(4)='temp';sv%varName(5)='sal'
    Globals%Constants%UseConstantWaterDensity=.true.
    Globals%Constants%UseConstantWaterDynamicViscosity=.false.
    Globals%Constants%MeanDensity=2._prec
    sv%state(:,4)=ieee_value(0._prec,ieee_quiet_nan)
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(all(abs(nu-0.001_prec/2._prec)<1e-18_prec),'computed NaN mu uses mean before division')
    call require(all(actualRho==2._prec),'invalid mu preserves chosen density')
    sv%state(:,4)=ieee_value(0._prec,ieee_positive_inf)
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(all(abs(nu-0.001_prec/2._prec)<1e-18_prec),'infinite input yields controlled mu')
    sv%state(:,4)=T
    ! Bad calculated density is corrected before dividing the independent mu by rho.
    Globals%Constants%UseConstantWaterDensity=.false.
    Globals%Constants%UseConstantWaterDynamicViscosity=.true.
    Globals%Constants%MeanDensity=1030._prec
    Globals%Constants%MeanKinematicViscosity=1.09e-6_prec
    Globals%Constants%MeanDynamicViscosity=0.001_prec
    sv%state(:,4)=[15._prec,ieee_value(0._prec,ieee_quiet_nan), &
                  ieee_value(0._prec,ieee_positive_inf),-1000._prec]
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(abs(actualRho(1)-baseRho(1))<1e-11_prec,'normal density remains calculated')
    call require(all(actualRho(2:)==1030._prec),'nonfinite or negative density uses mean')
    rho=water%seaWaterDensity(S,sv%state(:,4),Z)
    call require(maxval(abs(rho-actualRho))<1e-11_prec,'density function itself guards invalid calculated values')
    call require(maxval(abs(nu-0.001_prec/actualRho))<1e-16_prec,'nu uses the corrected density')
    call require(all(water%absoluteSeaWaterViscosity(S,T)==0.001_prec),'density guard leaves prescribed mu unchanged')
    sv%state(:,4)=1000._prec
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(all(actualRho==1030._prec),'large calculated density uses mean')
    call require(all(water%seaWaterDensity(S,sv%state(:,4),Z)==1030._prec), &
                 'density function itself guards large relative density')
    call require(all(abs(nu-0.001_prec/1030._prec)<1e-16_prec),'density correction precedes nu calculation')
    sv%state(:,4)=T
    Globals%Constants%MeanDensity=20000._prec
    call floater%DensityKViscosity_seawater(sv,nu,actualRho)
    call require(all(actualRho==20000._prec),'small relative calculated density uses mean')
    call require(all(water%seaWaterDensity(S,T,Z)==20000._prec), &
                 'density function itself guards small relative density')
    Globals%Constants%MeanDensity=1030._prec
    ! Actual buoyancy kernel: normal/large density contrast with and without T/S.
    allocate(buoySV%state(4,14),buoySV%varName(14))
    buoySV%varName(1)='x';buoySV%varName(2)='y';buoySV%varName(3)='z'
    buoySV%varName(4)='temp';buoySV%varName(5)='sal'
    buoySV%varName(6)='dist2bottom';buoySV%varName(7)='radius';buoySV%varName(8)='density'
    buoySV%varName(9)='volume';buoySV%varName(10)='area';buoySV%varName(11)='age'
    buoySV%varName(12)='rugosity';buoySV%varName(13)='dwz';buoySV%varName(14)='particulate'
    buoySV%state=0
    buoySV%state(:,3)=Z;buoySV%state(:,4)=T;buoySV%state(:,5)=S
    buoySV%state(:,7)=0.001_prec;buoySV%state(:,8)=[1000._prec,2000._prec,2200._prec,1100._prec]
    buoySV%state(:,9)=1e-8_prec;buoySV%state(:,10)=1e-6_prec;buoySV%state(:,11)=1._prec
    buoySV%state(:,12)=0.0025_prec;buoySV%state(:,13)=10._prec;buoySV%state(:,14)=1._prec
    Globals%Var%rugosityVar='rugosity';Globals%Var%dwz='dwz'
    Globals%Constants%Gravity%z=-9.81_prec
    Globals%Constants%MeanDensity=1100._prec
    Globals%Constants%MeanKinematicViscosity=1.09e-6_prec
    Globals%Constants%MeanDynamicViscosity=0.001_prec
    Globals%Constants%UseConstantWaterDynamicViscosity=.true.
    do i=0,2
        Globals%Constants%UseConstantWaterDensity=(i==1)
        if(i==2)then
            buoySV%varName(4)='unused1';buoySV%varName(5)='unused2'
        endif
        if(i==0)then
            usedRho=baseRho
        elseif(i==1)then
            usedRho=1100._prec
        else
            usedRho=ref_seaWaterDensity(0._prec*Z+36._prec,0._prec*Z+15._prec,Z)
        endif
        usedNu=0.001_prec/usedRho
        contrast=abs(1._prec-buoySV%state(:,8)/usedRho)
        if(i/=2)then
            where(contrast>=0.9_prec)
                contrast=abs(1._prec-buoySV%state(:,8)/1100._prec)
            endwhere
        endif
        diameter=(0.002_prec**3)*(9.81_prec*contrast/(usedNu**2))
        expectedSpeed=(usedNu/0.002_prec)*diameter* &
            (38.1_prec+0.93_prec*diameter**(4._prec/7._prec))**(-7._prec/8._prec)
        where(usedRho<buoySV%state(:,8))expectedSpeed=-expectedSpeed
        velocity=water%Buoyancy(buoySV,emptyBackground,0._prec)
        call require(maxval(abs(velocity(:,3)-expectedSpeed))<1e-12_prec, &
                     'buoyancy contrast guard matches expected fallback in each property mode')
    enddo
    ! Read explicit zeros through the real XML parser; they must not be discarded.
    path='water_properties.xml'
    open(unit=77,file='water_properties.xml',status='replace')
    write(77,'(a)') '<constants><UseConstantWaterDensity value="0"/>'// &
        '<UseConstantWaterDynamicViscosity value="0"/>'// &
        '<MeanDynamicViscosity value="0.002"/></constants>'
    close(77)
    call XMLReader%getFile(document,path)
    constants=>getDocumentElement(document);attribute='value'
    tag='UseConstantWaterDensity'
    call XMLReader%getNodeAttribute(constants,tag,attribute,text,found,.false.)
    call require(found,'XML parser retains explicit zero density flag')
    call Globals%Constants%setUseConstantWaterDensity(text)
    call require(.not.Globals%Constants%UseConstantWaterDensity,'XML zero disables density override')
    tag='UseConstantWaterDynamicViscosity'
    call XMLReader%getNodeAttribute(constants,tag,attribute,text,found,.false.)
    call require(found,'XML parser retains explicit zero viscosity flag')
    call Globals%Constants%setUseConstantWaterDynamicViscosity(text)
    call require(.not.Globals%Constants%UseConstantWaterDynamicViscosity,'XML zero disables viscosity override')
    tag='MeanDynamicViscosity'
    call XMLReader%getNodeAttribute(constants,tag,attribute,text,found,.false.)
    call require(found,'XML parser retains prescribed mu')
    call Globals%Constants%setMeanDynamicViscosity(text)
    call require(Globals%Constants%MeanDynamicViscosity==0.002_prec,'XML prescribed mu is read')
    tag='MissingFlag'
    call XMLReader%getNodeAttribute(constants,tag,attribute,text,found,.false.)
    call require(.not.found,'omitted XML flag is optional')
    call XMLReader%closeFile(document)
    print '(a,i0,a)', 'PASS: ',checks,' water-property checks'
contains
    subroutine require(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition)then
            print *, 'FAIL: ',message
            stop 1
        endif
        checks=checks+1
    end subroutine require
end program test_water_properties
