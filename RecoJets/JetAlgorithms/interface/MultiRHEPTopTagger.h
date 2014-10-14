//Adapted Version based on the work of Gregor.

#ifndef __MULTIR_TOPTAGGER_HH__
#define __MULTIR_TOPTAGGER_HH__

#include <vector>
#include <algorithm>  // for swap
#include <math.h>
#include "../interface/HEPTopTagger.h"

using namespace std;
using namespace fastjet;

// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
namespace external {


class MultiR_TopTagger {
public:
  MultiR_TopTagger(double max_fatjet_R,
		   double min_fatjet_R,
		   double step_R,
		   double multiR_threshold,
		   bool use_dR_max_triplet,
		   fastjet::ClusterSequence & cs, 
		   fastjet::PseudoJet & jet, 
		   double mtmass, double mwmass
		   );

  ~MultiR_TopTagger();

  //run tagger
  void run_tagger();

  const map<int,external::HEPTopTagger> & HTTagger() const {return _HEPTopTagger;}
  external::HEPTopTagger HTTagger(int i)  {return _HEPTopTagger[i];}
  
  const map<int,int> & n_small_fatjets() const {return _n_small_fatjets;}
  const double & Rmin() const {return _Rmin;}
  const double & mass_Rmin() const {return _mass_Rmin;}
  const double & pt_Rmin() const {return _pt_Rmin;}

  void set_max_subjet_mass(double x) {_subjet_mass = x;}
  void set_top_range(double top_range_min, double top_range_max) {_top_range[0] = top_range_min; _top_range[1] = top_range_max;}
  void set_f_W(double f_W) {_f_W = f_W;}
  void set_mass_ratio_cut(double mass_ratios_0, double mass_ratios_1, double mass_ratios_2) {_mass_ratios[0] = mass_ratios_0; _mass_ratios[1] = mass_ratios_1; _mass_ratios[2] = mass_ratios_2;}
  void set_nfilt(unsigned nfilt) {_n_filt = nfilt;}
  void set_Rfilt(double Rfilt) {_R_filt = Rfilt;}
  void set_debug(bool debug) {_debug = debug;}
 

private:
  const ClusterSequence * _cs;
  const PseudoJet *       _jet;
  double _mtmass, _mwmass;
  double _subjet_mass;
  map<int,external::HEPTopTagger> _HEPTopTagger;
  map<int,int> _n_small_fatjets;
  double _Rmin, _mass_Rmin, _pt_Rmin;
  double _mass_mean, _mass_width;
  double _top_range[2];
  unsigned _n_filt;
  double _R_filt;
  double _f_W;
  double _mass_ratios[3];
  double _max_fatjet_R, _min_fatjet_R, _step_R, _multiR_threshold;
  bool _use_dR_max_triplet;
  bool _debug;

  void UnclusterFatjets(const vector<fastjet::PseudoJet> & big_fatjets, vector<fastjet::PseudoJet> & small_fatjets, const ClusterSequence & cs, const double small_radius);

};
//--------------------------------------------------------------------
// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
};

#endif // __MULTIR_TOPTAGGER_HH__

