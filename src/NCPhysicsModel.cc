#include <iostream>
#include "NCPhysicsModel.hh"

#include "NCrystal/NCProcImpl.hh"
#include "NCrystal/NCException.hh"
#include "NCrystal/NCDefs.hh"

#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/NCSABScatter.hh"
#include "NCrystal/internal/NCSABFactory.hh"
#include "NCrystal/internal/NCSABExtender.hh"
#include "NCrystal/internal/NCSABIntegrator.hh"
#include "NCrystal/internal/NCSABUtils.hh"


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

NCP::PhysicsModel::PhysicsModel(const NC::Info& info)
  : m_proc(nullptr)
{
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto raw = info.getCustomSection( pluginNameUpperCase() );

  NC::ScatKnlData phononSab;
  phononSab.betaGridOptimised = true;

  
  phononSab.temperature = NC::Temperature{-1};
  // AtomMass should not have any impact to the result, as SABNullExtender 
  // is used for the energy range beyond the range of the single phonon sab 
  phononSab.elementMassAMU = NC::AtomMass{0.1}; 
  // The incoherent model in NCrystal is going to scale sab by ScatKnlData.SigmaBound.
  // However the bound scatting lengths are already included in the sab, so it should be unity. 
  phononSab.boundXS = NC::SigmaBound{1};


  NC::VectD *curField(nullptr);
  double num(0.);
  for(const auto& line : raw)
  {
    for(const auto& word : line)
    {
      if(NC::safe_str2dbl( word, num))
      {
        curField->push_back(num);
      }
      else
      {
        if(word=="alphagrid") {
          curField = &phononSab.alphaGrid;
        }
        else if(word=="betagrid") {
          curField = &phononSab.betaGrid;
        }
        else if(word=="sab_scaled") {
          curField = &phononSab.sab;
          phononSab.knltype = NC::ScatKnlData::KnlType::SCALED_SYM_SAB;
        }
        else if(word=="sab") {
          curField = &phononSab.sab;
          phononSab.knltype = NC::ScatKnlData::KnlType::SAB;
        }
        else if(word=="temperature")
        {
          double t(0);
          NC::safe_str2dbl(line[1], t);
          phononSab.temperature.set(t);
          std::cout << "temp " << t << std::endl;
          break;          
        }
        else
        {
          auto pos = word.find('r');
          if(pos==word.npos) {
            NCRYSTAL_THROW2(BadInput, "field name " << word << " is unknow ");  
          }  
          else
          {
            std::string digit = word.substr(0, pos);
            if(!NC::safe_str2dbl(digit, num))
              NCRYSTAL_THROW2(BadInput, "can not parse " << digit << " in " << word << " as double");    
            
            std::string times = word.substr(pos+1, word.size()-1);
            int rep(0);
            NC::safe_str2int(times, rep);
            for(unsigned i=0;i<rep;i++)
            {
              curField->push_back(num);
            }
          }
        }
      }
    }
  }

  // The code for calculating EmaxUpperBound is copied from the NC::validateScatKnlData function in NCScatKnlData.cc
  const double bmin = phononSab.betaGrid.front();
  const double amax = phononSab.alphaGrid.back();
  const double EmaxUpperBound = NC::constant_boltzmann*phononSab.temperature.get()*(bmin-amax)*(bmin-amax)/(4*amax);

  phononSab.suggestedEmax = NC::ncmin(EmaxUpperBound, 100.0);


  if(phononSab.alphaGrid.size()*phononSab.betaGrid.size()!=phononSab.sab.size())
    NCRYSTAL_THROW(BadInput, "data size error");

  NC::SABData sglphdata = NC::SABUtils::transformKernelToStdFormat(std::move(phononSab));

  auto sglphsab = NC::makeSO<NC::SABData>(std::move(sglphdata));
  auto sglphintegrator = NC::makeSO<NC::SAB::SABIntegrator>(std::move(sglphsab), nullptr, std::move(NC::makeSO<NC::SAB::SABNullExtender>()));
  auto sglphhelper = NC::makeSO<const NC::SAB::SABScatterHelper>(std::move(sglphintegrator->createScatterHelper()));
  auto sglphscatter = NC::makeSO<NC::SABScatter>(sglphhelper);
  NC::ProcImpl::ProcComposition::ComponentList components;
  components.push_back({1., sglphscatter});

  // if the line "auto sc_std = globalCreateScatter( cfg );"" in the NCPlugin Factory changes to
  // "auto sc_std = globalCreateScatter( cfg.modified("inelas=0") );", the lines followed show how we can create all the elastic
  // scattering processes from scratch. This method is potentially useful but may not be friendly for the final integration into the main repo. 
  // for ( auto& di : info.getDynamicInfoList() ) 
  // {
  //   auto di_vdos = dynamic_cast<const NC::DI_VDOS*>(di.get());
  //   if (di_vdos) 
  //   {
  //     double incohxs = di_vdos->atomData().incoherentXS().dbl();
  //     double scatteringxs = di_vdos->atomData().scatteringXS().dbl();
  //     double factor = incohxs/scatteringxs; 

  //     NC::ScatKnlData skd = NC::createScatteringKernel(di_vdos->vdosData(), 3,  0, NC::VDOSGn::TruncAndThinningChoices::Default,
  //                                 [factor](unsigned n) { 
  //                                   if(n>1) return 1.;
  //                                   else return factor; } );
      
  //     NC::SABData sabdata = NC::SABUtils::transformKernelToStdFormat(std::move(skd));
  //     components.push_back({di->fraction(), NC::makeSO<NC::SABScatter>(std::move(sabdata))});
      
  //   } 
  //   else
  //     NCRYSTAL_THROW(CalcError, "Not implemented error");
  // }

  NC::ProcImpl::ProcPtr procptr = NC::ProcImpl::ProcComposition::consumeAndCombine( std::move(components), NC::ProcessType::Scatter );

  auto rngproducer = NC::getDefaultRNGProducer();
  auto rng = rngproducer->produce();

  m_proc = std::make_shared<NC::Scatter>( std::move(rngproducer), std::move(rng), std::move(procptr));

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


