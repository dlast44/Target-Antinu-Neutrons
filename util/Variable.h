#ifndef VARIABLE_H
#define VARIABLE_H

//Includes from this package
#include "event/CVUniverse.h"
#include "util/SafeROOTName.h"
#include "util/Categorized.h"

//PlotUtils includes
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"

class Variable: public PlotUtils::VariableBase<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...)
    {
    }

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {

      std::map<int, std::string> BKGLabels = {{0, "NC"},
					       {1, "Wrong_Sign"}};

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Untested Changes                                                                                                                           //
      std::map<int, std::string> IntTypeLabels = {{1, "QE"},
						  {2, "RES"},
						  {3, "DIS"},
						  {8, "2p2h"}};

      //TODO: Get the actual ints correct
      std::map<int, std::string> TargetTypeLabels = {{1, "C"},
						     {2, "Fe"},
						     {3, "Pb"},
						     {8, "O"},
						     {9, "H"}};

      //TODO: Get the actual ints correct
      std::map<int, std::string> LeadBlobTypeLabels = {{1, "neut"},
						       {2, "mu"},
						       {3, "pi0"},
						       {8, "pim"},
						       {9, "pip"},
						       {10, "prot"},
						       {11, "None"}};
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      m_backgroundHists = new util::Categorized<Hist, int>((GetName() + "_background").c_str(),
							   GetName().c_str(), BKGLabels,
							   GetBinVec(), mc_error_bands);

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Untested Changes                                                                                                                           //
      m_SigIntTypeHists = new util::Categorized<Hist, int>((GetName() + "_sig").c_str(),
							   GetName().c_str(), IntTypeLabels,
							   GetBinVec(), mc_error_bands);

      m_SigTargetTypeHists = new util::Categorized<Hist, int>((GetName() + "_sig").c_str(),
							   GetName().c_str(), TargetTypeLabels,
							   GetBinVec(), mc_error_bands);

      m_SigLeadBlobTypeHists = new util::Categorized<Hist, int>((GetName() + "_sig").c_str(),
							   GetName().c_str(), LeadBlobTypeLabels,
							   GetBinVec(), mc_error_bands);

      m_BkgIntTypeHists = new util::Categorized<Hist, int>((GetName() + "_bkg").c_str(),
							   GetName().c_str(), IntTypeLabels,
							   GetBinVec(), mc_error_bands);

      m_BkgTargetTypeHists = new util::Categorized<Hist, int>((GetName() + "_bkg").c_str(),
							   GetName().c_str(), TargetTypeLabels,
							   GetBinVec(), mc_error_bands);

      m_BkgLeadBlobTypeHists = new util::Categorized<Hist, int>((GetName() + "_bkg").c_str(),
							   GetName().c_str(), LeadBlobTypeLabels,
							   GetBinVec(), mc_error_bands);
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      efficiencyNumerator = new Hist((GetName() + "_efficiency_numerator").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      efficiencyDenominator = new Hist((GetName() + "_efficiency_denominator").c_str(), GetName().c_str(), GetBinVec(), truth_error_bands);
      selectedSignalReco = new Hist((GetName() + "_selected_signal_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      selectedMCReco = new Hist((GetName() + "_selected_mc_reco").c_str(), GetName().c_str(), GetBinVec(), mc_error_bands);
      migration = new PlotUtils::Hist2DWrapper<CVUniverse>((GetName() + "_migration").c_str(), GetName().c_str(), GetBinVec(), GetBinVec(), mc_error_bands);
    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_backgroundHists;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Untested Changes                                                                                                                           //
    util::Categorized<Hist, int>* m_SigIntTypeHists;
    util::Categorized<Hist, int>* m_SigTargetTypeHists;
    util::Categorized<Hist, int>* m_SigLeadBlobTypeHists;
    util::Categorized<Hist, int>* m_BkgIntTypeHists;
    util::Categorized<Hist, int>* m_BkgTargetTypeHists;
    util::Categorized<Hist, int>* m_BkgLeadBlobTypeHists;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Hist* dataHist;
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;
    Hist* selectedSignalReco; //Effectively "true background subtracted" distribution for warping studies.
                              //Also useful for a bakground breakdown plot that you'd use to start background subtraction studies.
    Hist* selectedMCReco; //Treat the MC CV just like data for the closure test
    PlotUtils::Hist2DWrapper<CVUniverse>* migration;

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
      dataHist = new Hist((GetName() + "_data").c_str(), GetName().c_str(), GetBinVec(), data_error_bands);
    }

    void WriteData(TFile& file)
    {
      if (dataHist->hist) {
                dataHist->hist->SetDirectory(&file);
                dataHist->hist->Write();
      }
    }

    void WriteMC(TFile& file)
    {
      SyncCVHistos();
      file.cd();

      m_backgroundHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Untested Changes                                                                                                                           //
      m_SigIntTypeHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_SigTargetTypeHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_SigLeadBlobTypeHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_BkgIntTypeHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_BkgTargetTypeHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });

      m_BkgLeadBlobTypeHists->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?                                                                                           
                                    });
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      if(efficiencyNumerator)
      {
        efficiencyNumerator->hist->SetDirectory(&file); //TODO: Can I get around having to call SetDirectory() this many times somehow?
        efficiencyNumerator->hist->Write();
      }

      if(efficiencyDenominator)
      {
        efficiencyDenominator->hist->SetDirectory(&file);
        efficiencyDenominator->hist->Write();
      }

      if(migration)
      {
        migration->hist->SetDirectory(&file); 
        migration->hist->Write();
      }

      if(selectedSignalReco)
      {
        selectedSignalReco->hist->SetDirectory(&file);
        selectedSignalReco->hist->Write();
      }

      if(selectedMCReco)
      {
        selectedMCReco->hist->SetDirectory(&file);
        selectedMCReco->hist->Write((GetName() + "_data").c_str()); //Make this histogram look just like the data for closure tests
      }
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      m_backgroundHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //Untested Changes                                                                                                                           //
      m_SigIntTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_SigTargetTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_SigLeadBlobTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgIntTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgTargetTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      m_BkgLeadBlobTypeHists->visit([](Hist& categ) { categ.SyncCVHistos(); });
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
      if(selectedSignalReco) selectedSignalReco->SyncCVHistos();
      if(selectedMCReco) selectedMCReco->SyncCVHistos();
      if(migration) migration->SyncCVHistos();
    }
};

#endif //VARIABLE_H
