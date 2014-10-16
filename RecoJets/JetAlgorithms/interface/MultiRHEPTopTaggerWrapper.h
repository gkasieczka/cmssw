//----------------------------------------------------------------------
//  This file is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This file is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  The GNU General Public License is available at
//  http://www.gnu.org/licenses/gpl.html or you can write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------

#ifndef __MULTIRHEPTOPTAGGER_WRAPPER_HH__
#define __MULTIRHEPTOPTAGGER_WRAPPER_HH__

#include <fastjet/tools/TopTaggerBase.hh>
#include <fastjet/CompositeJetStructure.hh>
#include <sstream>

FASTJET_BEGIN_NAMESPACE

// Based on HEPTopTaggerWrapper.h

class MultiRHEPTopTaggerStructure;

class MultiRHEPTopTagger : public TopTaggerBase {
public:
  MultiRHEPTopTagger(double minSubjetPt, 
		     double minCandPt, 
		     double subjetMass, 
		     double muCut, 
		     int mode, 
		     double minCandMass, 
		     double maxCandMass, 
		     double massRatioWidth, 
		     double minM23Cut, 
		     double minM13Cut, 
		     double maxM13Cut,
		     double R_max,
		     double R_min) : minSubjetPt_(minSubjetPt),
    minCandPt_(minCandPt),
    subjetMass_(subjetMass),
    muCut_(muCut),
    mode_(mode),
    minCandMass_(minCandMass),
    maxCandMass_(maxCandMass),
    massRatioWidth_(massRatioWidth),
    minM23Cut_(minM23Cut),
    minM13Cut_(minM13Cut),
    maxM13Cut_(maxM13Cut),
    R_max_(R_max),
    R_min_(R_min)
  {}

  /// returns a textual description of the tagger
  virtual std::string description() const;

  /// runs the tagger on the given jet and
  /// returns the tagged PseudoJet if successful, or a PseudoJet==0 otherwise
  /// (standard access is through operator()).
  ///  \param jet   the PseudoJet to tag
  virtual PseudoJet result(const PseudoJet & jet) const;

  // the type of the associated structure
  typedef MultiRHEPTopTaggerStructure StructureType;

private: 
    double minSubjetPt_; // Minimal pT for subjets [GeV]
    double minCandPt_;   // Minimal pT to return a candidate [GeV]
 
    double subjetMass_; // Mass above which subjets are further unclustered
    double muCut_; // Mass drop threshold
    
    // HEPTopTagger Mode
    // 0: do 2d-plane, return candidate with delta m_top minimal
    // 1: return candidate with delta m_top minimal IF passes 2d plane
    // 2: do 2d-plane, return candidate with max dj_sum
    // 3: return candidate with max dj_sum IF passes 2d plane
    // 4: return candidate built from leading three subjets after unclustering IF passes 2d plane
    // Note: Original HTT was mode==1    
    int mode_; 

    // Top Quark mass window in GeV
    double minCandMass_;
    double maxCandMass_;
    
    double massRatioWidth_; // One sided width of the A-shaped window around m_W/m_top in %
    double minM23Cut_; // minimal value of m23/m123
    double minM13Cut_; // minimal value of atan(m13/m12)
    double maxM13Cut_; // maximal value of atan(m13/m12)

    double R_max_; // Maximal fatjet size to consider
    double R_min_; // Minimal fatjet size to consider
};


class MultiRHEPTopTaggerStructure : public CompositeJetStructure, public TopTaggerBaseStructure {
 public:
   /// ctor with pieces initialisation
   MultiRHEPTopTaggerStructure(const std::vector<PseudoJet>& pieces_in,
                  const JetDefinition::Recombiner *recombiner = 0) : CompositeJetStructure(pieces_in, recombiner),
    _top_mass(0.0),
    _unfiltered_mass(0.0),
    _pruned_mass(0.0),
    _fW(-1.),
    _mass_ratio_passed(-1),
    W_rec(recombiner), 
    rW_(){}
  
   // Return W subjet
   inline PseudoJet const & W() const{ 
     rW_ = join(_pieces[0], _pieces[1], *W_rec);
     return rW_;
   }
     
   // Return leading subjet in W
   inline PseudoJet  W1() const{
     assert(W().pieces().size()>0);
     return W().pieces()[0];
   }
       
   /// returns the second W subjet
   inline PseudoJet W2() const{
     assert(W().pieces().size()>1);
     return W().pieces()[1];
   }
 

   /// returns the non-W subjet
   /// It will have 1 or 2 pieces depending on whether the tagger has
   /// found 3 or 4 pieces
   inline const PseudoJet & non_W() const{ 
     return _pieces[2];
   }
 
   /// returns the candidate mass
   inline double top_mass() const {return _top_mass;}

   /// returns the unfiltered mass
   inline double unfiltered_mass() const {return _unfiltered_mass;}

   /// returns the pruned mass
   inline double pruned_mass() const {return _pruned_mass;}

   /// returns fW
   inline double fW() const {return _fW;}

   /// returns Rmin
   inline double R_min() const {return _Rmin;}

   /// returns if 2d-mass plane cuts were passed
   inline double mass_ratio_passed() const {return _mass_ratio_passed;}
    
 protected:
      double _top_mass;
      double _unfiltered_mass;
      double _pruned_mass;
      double _fW;      
      int _mass_ratio_passed;
      double _Rmin;
      
      const JetDefinition::Recombiner  * W_rec;
 
      mutable PseudoJet rW_;

   // allow the tagger to set these
   friend class HEPTopTagger;
   friend class MultiRHEPTopTagger;
 };


//------------------------------------------------------------------------
// description of the tagger
inline std::string MultiRHEPTopTagger::description() const{ 

  std::ostringstream oss;
  oss << "MultiRHEPTopTagger with: "
      << "minSubjetPt = " << minSubjetPt_ 
      << "minCandPt = " << minCandPt_ 
      << "subjetMass = " << subjetMass_ 
      << "muCut = " << muCut_ 
      << "mode = " << mode_ 
      << "minCandMass = " << minCandMass_ 
      << "maxCandMass = " << maxCandMass_ 
      << "massRatioWidth = " << massRatioWidth_ 
      << "minM23Cut = " << minM23Cut_ 
      << "minM13Cut = " << minM13Cut_
      << "Rmax = " << R_max_ 
      << "Rmin = " << R_min_ << std::endl;
  return oss.str();
}


FASTJET_END_NAMESPACE

#endif // __MULTIRHEPTOPTAGGER_HH__
