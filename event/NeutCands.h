//File: NeutCands.h
//Info: Neutron Candidate Class. Currently under development.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef NEUTCANDS_H
#define NEUTCANDS_H

#include "TVector3.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <bitset>

namespace NeutronCandidates{
  typedef std::map<std::string, std::vector<const char*>> intBranchMap;
  typedef std::map<std::string, std::vector<const char*>> doubleBranchMap;
  typedef std::map<std::string, std::vector<int>> intCandData;
  typedef std::map<std::string, std::vector<double>> doubleCandData;

  intBranchMap GetBranchIntMap();
  doubleBranchMap GetBranchDoubleMap();

  std::unordered_map<int,int> GetPDGBins();

  class NeutCand{
  private:
    //Currently only coding in the members that I actively use in MnvTgtNeutrons/particleCannon/nonMAT/interactiveMacros/Basic_Cuts_Try.cc
    
    int fID;
    int fIs3D;
    int fMCPID;
    int fTopMCPID;
    int fMCParentTrackID;
    int fMCParentPID;
    double fTotE;
    double fAngleToFP;
    TVector3 fEvtVtx;
    TVector3 fBegPos;
    TVector3 fEndPos;
    TVector3 fDirection;
    TVector3 fFlightPath;

    void init();

  public:
    //CTOR
    NeutCand();
    NeutCand(NeutronCandidates::intCandData candIntData, NeutronCandidates::doubleCandData candDoubleData, TVector3 vtx);

    //Here for use in Variable-like classes
    double GetDummyVar() const { return -999.; }

    int GetID() const { return fID; };
    int GetIs3D() const { return fIs3D; };
    int GetMCPID() const { return fMCPID; };
    int GetTopMCPID() const { return fTopMCPID; };
    int GetMCParentTrackID() const { return fMCParentTrackID; };
    int GetMCParentPID() const { return fMCParentPID; };
    
    double GetTotalE() const { return fTotE; };
    double GetAngleToFP() const { return fAngleToFP; };

    double GetPDGBin() const { return GetPDGBins()[fTopMCPID]; };
    double GetLength() const { return fDirection.Mag(); };
    double GetDEDX() const { 
      if (GetLength() > 0) return fTotE/GetLength(); 
      else return -1.0;
    };
    double GetVtxDist() const { return fFlightPath.Mag(); };
    double GetVtxZDist() const { return abs(fFlightPath.Z()); };
    double GetXPos() const { return fBegPos.X(); }
    double GetYPos() const { return fBegPos.Y(); }
    double GetZPos() const { return fBegPos.Z(); }

    TVector3 GetBegPos() const { return fBegPos; };
    TVector3 GetEndPos() const { return fEndPos; };
    TVector3 GetFlightPath() const { return fFlightPath; };
    TVector3 GetDirection() const { return fDirection; };
    TVector3 GetEvtVtx() const { return fEvtVtx; };
    std::bitset<4> GetClassifier();

    void SetID(std::vector<int> ID){ fID=ID.at(0); };
    void SetIs3D(std::vector<int> is3D){ fIs3D=is3D.at(0); };
    void SetMCPID(std::vector<int> MCPID){ fMCPID=MCPID.at(0); };
    void SetTopMCPID(std::vector<int> TopPID){ fTopMCPID=TopPID.at(0); };
    void SetMCParentTrackID(std::vector<int> ParentID){ fMCParentTrackID=ParentID.at(0); };
    void SetMCParentPID(std::vector<int> ParentPID){ fMCParentPID=ParentPID.at(0); };
    void SetTotalE(std::vector<double> TotE){ fTotE=TotE.at(0); };
    void SetEvtVtx(TVector3 EvtVtx){ fEvtVtx=EvtVtx; };
    //Move MULTI-LINE DEFINITIONS TO CPP...???
    void SetBegPos(std::vector<double> BegPos){
      fBegPos.SetXYZ(BegPos.at(0),BegPos.at(1),BegPos.at(2));
      fDirection = fEndPos-fBegPos;
      fFlightPath = fBegPos-fEvtVtx;
      if (fFlightPath.Mag() > 0 && fDirection.Mag() > 0){
	fAngleToFP = fFlightPath.Angle(fDirection);
      }
      else {
	fAngleToFP = -9999.0;
      }
    };
    void SetEndPos(std::vector<double> EndPos){
      fEndPos.SetXYZ(EndPos.at(0),EndPos.at(1),EndPos.at(2));
      fDirection = fEndPos-fBegPos;
      fAngleToFP = fFlightPath.Angle(fDirection);
      if (fFlightPath.Mag() > 0 && fDirection.Mag() > 0){
	fAngleToFP = fFlightPath.Angle(fDirection);
      }
      else {
	fAngleToFP = -9999.0;
      }
    };

    //DTOR
    virtual ~NeutCand() = default;
  };

  class NeutCands {
  private:
    int fNCands;
    int fIDmaxE;
    NeutCand fCandMaxE;
    std::map<int, NeutCand> fCands;

    void init();

  public:
    //CTORS
    NeutCands();
    NeutCands(std::vector<NeutCand> cands);
    NeutCands(std::map<int,intCandData> candsDataInt, std::map<int,doubleCandData> candDataDouble);

    //DTOR
    virtual ~NeutCands() = default;
    
    void SetCands(std::map<int, NeutCand> inCands){ 
      fCands=inCands; 
      fNCands=inCands.size();
    }

    int GetIDMaxE()  const { return fIDmaxE; };
    int GetNCands()  const { return fNCands; };
    NeutCand GetCandidate(int ID){ 
      if (fNCands == 0) return NeutCand();
      else return fCands[ID]; };
    NeutCand GetMaxCandidate() const { return fCandMaxE; };
    std::map<int, NeutCand> GetCandidates() const { return fCands; };
  };

}

class NeutronEvent{
 private:
  bool fIsSignal;
  int fIntType;
  int fTgtZ;
  NeutronCandidates::NeutCands fNeutCands;
 public:
  NeutronEvent() : fNeutCands() { }
  NeutronEvent(NeutronCandidates::NeutCands cands) { fNeutCands = cands; }

  bool IsSignal() const { return fIsSignal; }
  int GetIntType() const { return fIntType; }
  int GetTgtZ() const { return fTgtZ; }
  NeutronCandidates::NeutCand GetLeadingNeutCand() const { return fNeutCands.GetMaxCandidate(); }
  NeutronCandidates::NeutCands GetNeutCands() const { return fNeutCands; }

  void SetSignal(bool isSignal){ fIsSignal = isSignal; }
  void SetIntType(int intType){ fIntType = intType; }
  void SetTgtZ(int tgtZ){ fTgtZ = tgtZ; }

};

#endif
