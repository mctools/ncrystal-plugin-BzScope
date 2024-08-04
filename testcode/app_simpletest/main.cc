#include "NCPhysicsModel.hh"
#include <iostream>
#include "NCrystal/internal/NCMath.hh"//for NC::linspace

//Other includes (<random> is for the std::mt19937_64 example below):
#include <random>

namespace NC = NCrystal;

int main()
{


///////////////////////////////////////////////////////
  // Sanity check the NCrystal installation (optional) //
  ///////////////////////////////////////////////////////

  //This checks that the included NCrystal headers and the linked NCrystal
  //library are from the same release of NCrystal:
  NC::libClashDetect();

  ///////////////////////////////////////
  // Setup random generator (optional) //
  ///////////////////////////////////////

  //NCrystal already ships with a high quality generator, so this is done here
  //merely as an example for users who needs to use their own source of random
  //numbers. The example wraps the C++11 mt19937_64 generator and only
  //implements the actualGenerate function. Additional functions can be
  //overridden to enable more advanced capabilities concerning multi-threading,
  //RNG state manipulation, etc. See the NCRNG.hh file for details. Of course,
  //if we had NOT registered a custom RNG source here, NCrystal would be using
  //the built-in source, which IS multi-thread safe (due to usage of jump-ahead
  //features).

  class CustomRNG : public NC::RNGStream {
    std::mt19937_64 m_gen;
  protected:
    double actualGenerate() override { return NC::randUInt64ToFP01(static_cast<uint64_t>(m_gen())); }
    //For the sake of example, we wrongly claim that this generator is safe and
    //sensible to use multithreaded (see NCRNG.hh for how to correctly deal with
    //MT safety, RNG states, etc.):
    bool useInAllThreads() const override { return true; }
  };

  //The NCrystal makeSO function is similar to std::make_shared
  //and should be used instead of raw calls to new and delete:
  auto rng = NC::makeSO<CustomRNG>();

  //Register:
  NC::setDefaultRNG(rng);

  //////////////////////////////////////
  // Create and use aluminium powder: //
  //////////////////////////////////////

  auto pc = NC::createScatter( "/home/caixx/git/ncplugin-BzScope/testcode/data/custom.ncmat" );



  return 0;
}
