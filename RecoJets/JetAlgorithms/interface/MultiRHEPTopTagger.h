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
		   const fastjet::ClusterSequence & cs, 
		   const fastjet::PseudoJet & jet, 
		   double mtmass, double mwmass
		   );

  ~MultiR_TopTagger();

  //run tagger
  void run_tagger();

  // Return the candidate (and some properties) at R=R_min
  HEPTopTagger cand_Rmin(){return _HEPTopTagger[_Rmin];}
  const int & Rmin_raw() const {return _Rmin;}
  const double Rmin() const {return _Rmin/10.;}
  const double & mass_Rmin() const {return _mass_Rmin;}
  const double & pt_Rmin() const {return _pt_Rmin;}
  
  double R_min_exp(double x) {return _r_min_exp_function(x);}

  // Access to all candidates and number-of-small-fatjets
  HEPTopTagger HTTagger(int i)  {return _HEPTopTagger[i];}
  const double n_small_fatjets(int i) {return _n_small_fatjets[i];}

  void set_mode(int mode) {_mode = mode;}

  void set_max_subjet_mass(double x) {_subjet_mass = x;}
  void set_mass_drop_threshold(double x) {_mass_drop_threshold = x;}
  void set_minpt_subjet(double x) {_minpt_subjet = x;}
  void set_minpt_tag(double x) {_minpt_tag = x;}
  void set_top_range(double top_range_min, double top_range_max) {_top_range[0] = top_range_min; _top_range[1] = top_range_max;}
  void set_f_W(double f_W) {_f_W = f_W;}
  void set_mass_ratio_cut(double mass_ratios_0, double mass_ratios_1, double mass_ratios_2) {_mass_ratios[0] = mass_ratios_0; _mass_ratios[1] = mass_ratios_1; _mass_ratios[2] = mass_ratios_2;}
  void set_nfilt(unsigned nfilt) {_n_filt = nfilt;}
  void set_Rfilt(double Rfilt) {_R_filt = Rfilt;}
  void set_debug(bool debug) {_debug = debug;}
  void set_r_min_exp_function(double (*f)(double)) {_r_min_exp_function = f;}
 

private:
  const ClusterSequence * _cs;
  const PseudoJet *       _jet;
  double _mtmass, _mwmass;
  double _mass_drop_threshold;
  double _subjet_mass;
  double _minpt_tag;
  double _minpt_subjet;
  map<int,external::HEPTopTagger> _HEPTopTagger;
  map<int,int> _n_small_fatjets;
  int _Rmin;
  int _mode;
  double _mass_Rmin, _pt_Rmin;
  double _mass_mean, _mass_width;
  double _top_range[2];
  unsigned _n_filt;
  double _R_filt;
  double _f_W;
  double _mass_ratios[3];
  double _max_fatjet_R, _min_fatjet_R, _step_R, _multiR_threshold;
  bool _use_dR_max_triplet;
  bool _debug;
  double (*_r_min_exp_function)(double);

  void UnclusterFatjets(const vector<fastjet::PseudoJet> & big_fatjets, vector<fastjet::PseudoJet> & small_fatjets, const ClusterSequence & cs, const double small_radius);

};
//--------------------------------------------------------------------
// Do not change next line, it's needed by the sed-code that makes the tagger CMSSW-compatible.
};

#endif // __MULTIR_TOPTAGGER_HH__

