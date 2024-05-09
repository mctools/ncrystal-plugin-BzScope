#include <iostream>
#include "NCPhysicsModel.hh"

#include "NCrystal/NCProcImpl.hh"
#include "NCrystal/NCException.hh"

#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/NCSABScatter.hh"
#include "NCrystal/internal/NCSABFactory.hh"
#include "NCrystal/internal/NCSABExtender.hh"
#include "NCrystal/internal/NCSABIntegrator.hh"


bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  if(!info.hasDynamicInfo())
    return false;

  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;

}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
{
  return PhysicsModel(info);
}

NCP::PhysicsModel::PhysicsModel( const NC::Info& info)
  : m_proc(nullptr)
{
    //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <sigmavalue> <wavelength threshold value>
  //

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Verify we have exactly one line and two words:
  if ( data.size() != 1 || data.at(0).size()!=2 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be two numbers on a single line");

  //Parse and validate values:
  double sigma, lambda_cutoff;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), sigma )
       || ! NC::safe_str2dbl( data.at(0).at(1), lambda_cutoff )
       || ! (sigma>0.0) || !(lambda_cutoff>=0.0) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" section (should be two positive floating point values)" );

  

  NC::ProcImpl::ProcComposition::ComponentList components;

  for ( auto& di : info.getDynamicInfoList() ) 
  {
    auto di_vdos = dynamic_cast<const NC::DI_VDOS*>(di.get());
    if (di_vdos) 
    {
      double incohxs = di_vdos->atomData().incoherentXS().dbl();
      double scatteringxs = di_vdos->atomData().scatteringXS().dbl();
      double factor = incohxs/scatteringxs;

      NC::ScatKnlData skd = NC::createScatteringKernel(di_vdos->vdosData(), 3,  0, NC::VDOSGn::TruncAndThinningChoices::Default); 
                                  // [factor](unsigned n) { 
                                  //   // std::cout << "order " << n  << "\n";
                                  //   if(n>1) return 1.;
                                  //   else return factor; } );

      NC::SABData sabdata(std::move(skd.alphaGrid), std::move(skd.betaGrid), std::move(skd.sab), skd.temperature, 
      skd.boundXS, skd.elementMassAMU);
      // components.push_back({di->fraction(), NC::makeSO<NC::SABScatter>(std::move(sabdata))});


      // or
      auto shrdsab = NC::makeSO<NC::SABData>(std::move(sabdata));
      auto integrator = NC::makeSO<NC::SAB::SABIntegrator>(std::move(shrdsab), nullptr, std::move(NC::makeSO<NC::SAB::SABNullExtender>()));
      auto helper = NC::makeSO<const NC::SAB::SABScatterHelper>(std::move(integrator->createScatterHelper()));
      auto scatter = NC::makeSO<NC::SABScatter>(helper);
      components.push_back({di->fraction(), (scatter)});

    } 
    else
      NCRYSTAL_THROW(CalcError, "Not implemented error");
  }

  // NC::SAB::SABScatterHelper a = NC::SAB::createScatterHelper( shared_obj<const NC::SABData> data,
  //                                                             std::shared_ptr<const VectD> energyGrid );
  // NC::SABXSProvider &pro = a.xsprovider;
  // auto en = pro.internalEGrid();
  // auto xs = pro.internalXSGrid();
  // pro.set(en, xs, NC::SABNullExtender());


  // std::unique_ptr<const NC::SAB::SABScatterHelper> NC::SAB::createScatterHelper( shared_obj<const NC::SABData> data,
  //                                                                              std::shared_ptr<const VectD> energyGrid )
  // {
  //   nc_assert(!!data);
  //   SABIntegrator si(data,energyGrid.get());
  //   auto sh = si.createScatterHelper();
  //   return std::make_unique<SABScatterHelper>(std::move(sh));
  // }



  // low order phonon scattering 

  // For the low order coherent scattering, the constructor SABScatter( SAB::SABScatterHelper&& ) should be used to control the extender,
  // which in this case should always return zero.
    
  // SABScatter( SAB::SABScatterHelper&& );

  // SABScatterHelper( SABXSProvider&& xp,
  //             SABSampler&&sp,
  //             Optional<std::string> json = NullOpt )

  // SABXSProvider( VectD&& egrid, VectD&& xsvals,
  //           std::shared_ptr<const SAB::SABExtender> );
  // where SABExtender should be SABNullExtender() and
  //     SABSampler( Temperature temperature,
  //                 VectD&& egrid,
  //                 std::vector<std::unique_ptr<SABSamplerAtE>>&&,
  //                 std::shared_ptr<const SAB::SABExtender>,
  //                 double xsAtEmax,
  //                 EGridMargin );





  NC::ProcImpl::ProcPtr procptr = NC::ProcImpl::ProcComposition::consumeAndCombine( std::move(components), NC::ProcessType::Scatter );

  auto rngproducer = NC::getDefaultRNGProducer();
  auto rng = rngproducer->produce();

  m_proc = std::make_shared<NC::Scatter>( std::move(rngproducer), std::move(rng), std::move(procptr));


  //Important note to developers who are using the infrastructure in the
  //testcode/ subdirectory: If you change the number or types of the arguments
  //for the constructor here, you should make sure to perform a corresponding
  //change in three files in the testcode/ directory: _cbindings.py,
  //__init__.py, and NCForPython.cc - that way you can still instantiate your
  //model directly from your python test code).


  // /////////////////////////////////////////////////////////////////////////////
  // // Expand a vibrational density of state (VDOS) spectrum into a full-blown //
  // // scattering kernel.                                                      //
  // The returned type is ScatKnlData::KnlType::SAB;
  // /////////////////////////////////////////////////////////////////////////////


  // using ScaleGnContributionFct = std::function<double(unsigned)>;
  // ScatKnlData createScatteringKernel( const VDOSData&,
  //                                     unsigned vdosluxlvl = 3,//0 to 5, affects binning, Emax, etc.
  //                                     double targetEmax = 0.0,//if 0, will depend on luxlvl. Error if set to unachievable value.
  //                                     const VDOSGn::TruncAndThinningParams ttpars = VDOSGn::TruncAndThinningChoices::Default,
  //                                     ScaleGnContributionFct = nullptr,
  //                                     Optional<unsigned> override_max_order = NullOpt );

  // //The ScaleGnContributionFct argument can be used to apply a scale factor to
  // //the Gn contribution, at the point where it is used to add a contribution
  // //into S(alpha,beta). This can for instance be used in the scenario where the
  // //coherent single-phonon contribution is modelled elsewhere, and therefore the
  // //G1 contribution should be reduced to take out sigma_coherent. In this case,
  // //the scaling function should return sigma_incoh/(sigma_coh+sigma_incoh) for
  // //the argument n==1, and 1.0 for n>1.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //2. need to convert ScatKnlData to DI_ScatKnl or SABData

  //   struct ScatKnlData : private MoveOnly {
  //   VectD alphaGrid, betaGrid, sab;
  //   Temperature temperature = Temperature{-1.0};
  //   SigmaBound boundXS = SigmaBound{-1.0};
  //   AtomMass elementMassAMU = AtomMass{-1.0};
  //   enum class KnlType { SAB,        //Standard S(alpha,beta), fields have same meaning as on SABData.
  //                        SCALED_SAB, //Values in sab table are actually S'(alpha,beta)=S(alpha,beta)*exp(beta/2)
  //                        SCALED_SYM_SAB,//Same + S'(alpha,beta) is an even function in beta, and the kernel is
  //                        //only specified directly for non-negative beta values.
  //                        SQW,//alpha/beta/sab values are actually Q/omega/S(Q,omega) values.
  //                        Unspecified };
  //   KnlType knltype = KnlType::Unspecified;
  //   double suggestedEmax = 0;//optional, 0 means not available.
  //   bool betaGridOptimised = false;//set to true if beta grid is known to not require thickening before sampling.
  // };


    //Construct from SABData and (optionally) an energy grid. The first
    //constructor takes a DI_ScatKnl object from the DynamicInfo list on an an
    //Info object, along with (optionally) the unique id of the Info object. If
    //an energy grid is not supplied, a reasonable default value will be
    //used. When possible, caching will be enabled by default - making sure that
    //multiple SABScatter instances based on the same input object will avoid
    //duplicated resource consumption.
    //
    //The vdoslux parameter has no effect if input is not a VDOS. The same goes
    //for the special vdos2sabExcludeFlag parameter (the meaning of which is
    //documented in NCDynInfoUtils.hh).


    //     //Constructors etc. (all expensive operations forbidden):
      // SABData( VectD&& alphaGrid, VectD&& betaGrid, VectD&& sab,
      //          Temperature temperature, SigmaBound boundXS, AtomMass elementMassAMU,
      //          double suggestedEmax = 0 );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 3. create the model

  // DI_ScatKnl is the base of DI_VDOS

  // SABScatter( const DI_ScatKnl&,
  //             unsigned vdoslux = 3,
  //             bool useCache = true,
  //             uint32_t vdos2sabExcludeFlag = 0 );
  // SABScatter( SABData &&,
  //             const VectD& energyGrid = VectD() );
  // SABScatter( shared_obj<const SABData>,
  //             std::shared_ptr<const VectD> energyGrid = nullptr );
  // explicit SABScatter( shared_obj<const SAB::SABScatterHelper> );
  // explicit SABScatter( std::unique_ptr<const SAB::SABScatterHelper> );
  // SABScatter( SAB::SABScatterHelper&& );



}

double NCP::PhysicsModel::calcCrossSection( double neutron_ekin ) const
{
  return m_proc->crossSectionIsotropic(NC::NeutronEnergy{neutron_ekin}).dbl();
}

NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const
{
  ScatEvent result;
  auto res = m_proc->sampleScatterIsotropic(NC::NeutronEnergy{neutron_ekin});
  result.ekin_final = res.ekin.dbl();
  result.mu = res.mu.dbl();
  return result;
}


