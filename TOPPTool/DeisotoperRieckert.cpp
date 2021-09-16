// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lucas Rieckert $
// $Authors: Lucas Rieckert $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <string>
#include <iterator>
#include <sstream>
#include <numeric>
#include <OpenMS/FILTERING/DATAREDUCTION/DeisotoperRieckert.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <time.h>
#include <chrono>

using namespace std;

namespace OpenMS
{
  ///Funktion zum düchführen der verschiedenen Algorithmentests, wie true positiv ect. Berechnung und die Überprüfung wie gut die Algorithmen mit manuel gewählten Isotopendistributionen zurecht kamen
  ///Da die ersten Test keine guten Ergebnisse lieferten passte ich die Funktion für den letzten Test an, der auskommentierte Teil diente den ersten Test und dient nur zur Veranschaulichung des Vorgehens und kann ohne Anpassungen nicht mehr ausgeführt werden.
  ///Der Test der manuellen Daten ist der nicht auskommentierte Teil. Dieser ist funktional.
  void DeisotoperRieckert::dataAnalysis(std::string mzml, std::string featureXMLin)
  {

    DTA2DFile file_builder;

    String inputMZMl = mzml;//"/buffer/ag_bsc/pmsb_2021/Rieckert/data/centroided-sample-human-orbi2.mzML";
    String raw_MZml = "/buffer/ag_bsc/pmsb_2021/Rieckert/data/raw-sample-human-orbi2.mzML";
    PeakMap exp_deiso;
    PeakMap exp_data;
    PeakMap exp_data_raw;
    PeakMap exp_resampled;
    PeakMap exp_extra;
    PeakMap exp_extra_tmp;

    PeakMap hardklor;
    PeakMap hardklor3;
    PeakMap decon;
    PeakMap openDeiso;

    MSExperiment hardklor_true_pos;
    MSExperiment hardklor_false_pos;
    MSExperiment hardklor_true_neg;
    MSExperiment hardklor_false_neg;

    MSExperiment hardklor3_true_pos;
    MSExperiment hardklor3_false_pos;
    MSExperiment hardklor3_true_neg;
    MSExperiment hardklor3_false_neg;

    MSExperiment decon_true_pos;
    MSExperiment decon_false_pos;
    MSExperiment decon_true_neg;
    MSExperiment decon_false_neg;

    MSExperiment openDeiso_true_pos;
    MSExperiment openDeiso_false_pos;
    MSExperiment openDeiso_true_neg;
    MSExperiment openDeiso_false_neg;

    bool peak_found;
    double mz;
    double rt;

    int all_peaks = 0;
    int all_peaks_raw = 0;
    MSExperiment non_mono_peaks;
    MSExperiment mono_peaks;

    MSExperiment non_mono_peaks_raw;
    MSExperiment mono_peaks_raw;

    MSSpectrum buffer_found;
    MSSpectrum buffer_not_found;

    MzMLFile mzFile;

    ///Datenimporte für die ersten Tests
   /*
    mzFile.load(inputMZMl, exp_deiso);
    mzFile.load(inputMZMl, exp_data);
    mzFile.load(raw_MZml,exp_data_raw);
    mzFile.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/resampled-raw-sample-human-orbi2.mzML", exp_resampled);
    mzFile.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/P4 dependent replicate forward 2.mzML", exp_extra);
    mzFile.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/MC_bud1_Chymotrypsin_1.mzML", exp_extra);
    file_builder.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/old2/filtPeaksRT.dta2d", hardklor);
    file_builder.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filtPeaksRT.dta2d", hardklor3);
    file_builder.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/Decon2LS.dta2d", decon);
    FeatureXMLFile featureXML;
    FeatureMap featureMap;
    featureXML.load(featureXMLin, featureMap); //"/buffer/ag_bsc/pmsb_2021/Rieckert/data/featureMap-sample-human-orbi2.featureXML"
    */


    ///Test mit den manuel gewählten Isotopendistributionen.

    ///Laden des getrimmten Mzml files für den opendeisotoper
    mzFile.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/trimmed/MC_bud3_Trypsin_1_trimmed_1412_92.mzML", exp_extra_tmp);

    PeakPickerHiRes pp;
    pp.pickExperiment(exp_extra_tmp,exp_extra);
    mzFile.store("/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/trimmed/MC_bud3_Trypsin_1_picked.mzML",exp_extra);

    //mzFile.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/trimmed/MC_bud3_Trypsin_1_picked.mzML",exp_extra);

    vector<MSSpectrum> spectra = exp_extra.getSpectra();

    double fragment_tolerance = 5.0;
    bool fragment_unit_ppm = true;
    int min_charge = 1;
    int max_charge = 10;
    bool keep_only_deisotoped = true;
    unsigned int min_isopeaks = 2;
    unsigned int max_isopeaks = 5;
    bool make_single_charged = false;
    bool annotate_charge = true;
    bool annotate_iso_peak_count = true;
    bool use_decreasing_model = true;
    unsigned int start_intensity_check = 1;
    bool add_up_intensity = false;

    exp_extra.updateRanges();

    ///OpenMS Deisotoper Run

    for (int i = 0; i < spectra.size(); i++)
    {

      Deisotoper::deisotopeAndSingleCharge(spectra[i],
                                           fragment_tolerance,
                                           fragment_unit_ppm,
                                           min_charge,
                                           max_charge,
                                           keep_only_deisotoped,
                                           min_isopeaks,
                                           max_isopeaks,
                                           make_single_charged,
                                           annotate_charge,
                                           annotate_iso_peak_count,
                                           use_decreasing_model,
                                           start_intensity_check,
                                           add_up_intensity);

    }

    openDeiso.setSpectra(spectra);
    openDeiso.updateRanges();


    vector<double> rt_buffer;
    vector<double> hardklor_peaks;
    vector<double> decon_peaks;
    vector<double> mono_gut;
    vector<double> mono_overlap;
    vector<double> hardklor_found_gut;
    vector<double> hardklor_found_overlap;
    vector<double> hardklor_not_found_gut;
    vector<double> hardklor_not_found_overlap;
    vector<double> decon_found_gut;
    vector<double> decon_found_overlap;
    vector<double> decon_not_found_gut;
    vector<double> decon_not_found_overlap;
    vector<double> open_found_gut;
    vector<double> open_found_overlap;
    vector<double> open_not_found_gut;
    vector<double> open_not_found_overlap;

    int hardklor_hits_gut=0;
    int hardklor_hits_overlap=0;
    int hardklor_no_hits_gut=0;
    int hardklor_no_hits_overlap=0;
    int decon_hits_gut=0;
    int decon_hits_overlap=0;
    int decon_no_hits_gut=0;
    int decon_no_hits_overlap=0;
    int open_hits_gut=0;
    int open_hits_overlap=0;
    int open_no_hits_gut=0;
    int open_no_hits_overlap=0;

    bool found;


    mono_gut.push_back(492.79514692);
    mono_gut.push_back(495.26890796);
    mono_gut.push_back(514.96538407);
    mono_gut.push_back(522.9216);
    mono_gut.push_back(532.5962);
    mono_gut.push_back(537.2940);
    mono_gut.push_back(539.2633);
    mono_gut.push_back(540.9581);
    mono_gut.push_back(543.2627);
    mono_gut.push_back(547.9644);
    mono_gut.push_back(551.7821);
    mono_gut.push_back(558.7944);
    mono_gut.push_back(666.8547);
    mono_gut.push_back(686.3636);
    mono_gut.push_back(690.8865);
    mono_gut.push_back(692.8922);
    mono_gut.push_back(694.8587);
    mono_gut.push_back(703.8440);
    mono_gut.push_back(708.3385);
    mono_gut.push_back(718.3331);
    mono_gut.push_back(726.3143);
    mono_gut.push_back(744.7236);
    mono_gut.push_back(761.4480);
    mono_gut.push_back(771.9433);
    mono_gut.push_back(775.3913);

    mono_overlap.push_back(502.18333738);
    mono_overlap.push_back(502.26224180);
    mono_overlap.push_back(503.78034266);
    mono_overlap.push_back(504.28447026);
    mono_overlap.push_back(505.61820004);
    mono_overlap.push_back(505.95209791);
    mono_overlap.push_back(507.96968839);
    mono_overlap.push_back(508.27628473);
    mono_overlap.push_back(509.26567794);
    mono_overlap.push_back(509.30255974);
    mono_overlap.push_back(517.7848);
    mono_overlap.push_back(518.9310);
    mono_overlap.push_back(639.3406);
    mono_overlap.push_back(639.6167);
    mono_overlap.push_back(711.8659);
    mono_overlap.push_back(712.3370);
    mono_overlap.push_back(733.8353);
    mono_overlap.push_back(733.9000);
    mono_overlap.push_back(821.4455);
    mono_overlap.push_back(823.3581);

    ifstream infileHardklor("/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/handpick/foundPeaksHardklor_1412_92_picked_c10_none_v2_d5.txt");
    ifstream infileDecon("/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/handpick/foundPeaksDecon_1412_92.txt");

    for(string line = ""; getline(infileHardklor,line); )
    {
      hardklor_peaks.push_back(stod(line));
    }

    for(string line = ""; getline(infileDecon, line); )
    {
      decon_peaks.push_back(stod(line));
    }
    infileDecon.close();
    infileHardklor.close();

    cout << "Tool output gelesen" <<endl;

    for(int i = 0; i < hardklor_peaks.size(); ++i)
    {
      for(int j = 0; j < mono_gut.size(); ++j)
      {
        if(abs(hardklor_peaks[i] - mono_gut[j])<0.01)
        {
          hardklor_hits_gut+=1;
          hardklor_found_gut.push_back(mono_gut[j]);
        }
      }

      for(int k = 0; k < mono_overlap.size(); ++k)
      {
        if(abs(hardklor_peaks[i] - mono_overlap[k])<0.01)
        {
          hardklor_hits_overlap+=1;
          hardklor_found_overlap.push_back(mono_overlap[k]);
        }
      }
    }

    for(int i =0; i< mono_gut.size(); ++i)
    {
      found=false;
      for(int j = 0; j < hardklor_found_gut.size(); ++j)
      {
        if(abs(hardklor_found_gut[j] - mono_gut[i]) <0.01)
        {
          found = true;
        }
      }
      if(!found)
      {
        hardklor_not_found_gut.push_back(mono_gut[i]);
        hardklor_no_hits_gut+=1;
      }
    }

    for(int i =0; i< mono_overlap.size(); ++i)
    {
      found=false;
      for(int j = 0; j < hardklor_found_overlap.size(); ++j)
      {
        if(abs(hardklor_found_overlap[j] - mono_overlap[i]) <0.01)
        {
          found = true;
        }
      }
      if(!found)
      {
        hardklor_not_found_overlap.push_back(mono_overlap[i]);
        hardklor_no_hits_overlap+=1;
      }
    }

    cout << "Hardklor done!" << endl;

    for(int i = 0; i<decon_peaks.size(); ++i)
    {
      cout<<decon_peaks[i] << " , " ;
    }
    cout << endl;

    for(int i = 0; i < decon_peaks.size(); ++i)
    {
      for(int j = 0; j < mono_gut.size(); ++j)
      {
        if(abs(decon_peaks[i] - mono_gut[j])<0.01)
        {
          cout << "Hit bei DP: " << decon_peaks[i] << " MONOGUT : " << mono_gut[j]<< endl;
          decon_hits_gut+=1;
          decon_found_gut.push_back(mono_gut[j]);
        }
      }

      for(int k = 0; k < mono_overlap.size(); ++k)
      {
        if(abs(decon_peaks[i] - mono_overlap[k])<0.01)
        {
          cout << "Hit bei DP: " << decon_peaks[i] << " MONOOVERLAP: " << mono_overlap[k]<< endl;
          decon_hits_overlap+=1;
          decon_found_overlap.push_back(mono_overlap[k]);
        }
      }
    }

    for(int i =0; i< mono_gut.size(); ++i)
    {
      found=false;
      for(int j = 0; j < decon_found_gut.size(); ++j)
      {
        if(abs(decon_found_gut[j] - mono_gut[i]) <0.01)
        {
          found = true;
        }
      }
      if(!found)
      {
        decon_not_found_gut.push_back(mono_gut[i]);
        decon_no_hits_gut+=1;
      }
    }

    for(int i =0; i< mono_overlap.size(); ++i)
    {
      found=false;
      for(int j = 0; j < decon_found_overlap.size(); ++j)
      {
        if(abs(decon_found_overlap[j] - mono_overlap[i]) <0.01)
        {
          found = true;
        }
      }
      if(!found)
      {
        decon_not_found_overlap.push_back(mono_overlap[i]);
        decon_no_hits_overlap+=1;
      }
    }

    cout << "decon done!" << endl;


    for(auto spec:openDeiso)
    {
      for(auto peak:spec)
      {
        //cout << peak.getMZ()<<": ";
        for(int j = 0; j < mono_gut.size(); j++)
        {
          if(abs(peak.getMZ() - mono_gut[j])<0.01)
          {
            open_hits_gut+=1;
            open_found_gut.push_back(mono_gut[j]);
            //cout<< peak.getMZ()<<": "<<j<<"!!!!,";
            break;
          }
        }

        for(int k = 0; k < mono_overlap.size(); k++)
        {
          if(abs(peak.getMZ() - mono_overlap[k])<0.01)
          {
            open_hits_overlap+=1;
            open_found_overlap.push_back(mono_overlap[k]);
            //cout<< peak.getMZ()<<": "<<k<<"!!!!!,";
            break;
          }
        }
      }
    }

    for(int i =0; i< mono_gut.size(); ++i)
    {
      found=false;
      for(int j = 0; j < open_found_gut.size(); ++j)
      {
        if(abs(open_found_gut[j] - mono_gut[i]) <0.01)
        {
          found = true;
        }
      }
      if(!found)
      {
        open_not_found_gut.push_back(mono_gut[i]);
        open_no_hits_gut+=1;
      }
    }

    for(int i =0; i< mono_overlap.size(); ++i)
    {
      found=false;
      for(int j = 0; j < open_found_overlap.size(); ++j)
      {
        if(abs(open_found_overlap[j] - mono_overlap[i]) <0.01)
        {
          found = true;
        }
      }
      if(!found)
      {
        open_not_found_overlap.push_back(mono_overlap[i]);
        open_no_hits_overlap+=1;
      }
    }

    cout << "open done!" << endl;

    cout << "Hardklor hat von " << mono_gut.size() << " guten Isotopenmuster " << hardklor_hits_gut << " erfolgreich gefunden." <<endl;
    cout << "Hardklor hat von " << mono_overlap.size() << " ueberlappenden Isotopenmuster " << hardklor_hits_overlap << " erfolgreich gefunden." <<endl;

    cout << "Die Guten sind: ";
    for(int i = 0; i< hardklor_found_gut.size();++i)
    {
      cout << hardklor_found_gut[i]<<", ";
    }
    cout<< endl;

    cout << "Die ueberlappenden sind: ";
    for(int i = 0; i< hardklor_found_overlap.size();++i)
    {
      cout << hardklor_found_overlap[i]<<", ";
    }
    cout<< endl;

    cout << "Hardklor hat von " << mono_gut.size() << " guten Isotopenmuster " << hardklor_no_hits_gut << " nicht gefunden." <<endl;
    cout << "Hardklor hat von " << mono_overlap.size() << " ueberlappenden Isotopenmuster " << hardklor_no_hits_overlap << " nicht gefunden." <<endl;

    cout << "Die Guten sind: ";
    for(int i = 0; i< hardklor_not_found_gut.size();++i)
    {
      cout << hardklor_not_found_gut[i]<<", ";
    }
    cout<< endl;

    cout << "Die ueberlappenden sind: ";
    for(int i = 0; i< hardklor_not_found_overlap.size();++i)
    {
      cout << hardklor_not_found_overlap[i]<<", ";
    }
    cout<< endl;


    cout << "Decon2LS hat von " << mono_gut.size() << " guten Isotopenmuster " << decon_hits_gut << " erfolgreich gefunden." <<endl;
    cout << "Decon2LS hat von " << mono_overlap.size() << " ueberlappenden Isotopenmuster " << decon_hits_overlap << " erfolgreich gefunden." <<endl;

    cout << "Die Guten sind: ";
    for(int i = 0; i< decon_found_gut.size();++i)
    {
      cout << decon_found_gut[i]<<", ";
    }
    cout<< endl;

    cout << "Die ueberlappenden sind: ";
    for(int i = 0; i< decon_found_overlap.size();++i)
    {
      cout << decon_found_overlap[i]<<", ";
    }
    cout<< endl;

    cout << "Decon2LS hat von " << mono_gut.size() << " guten Isotopenmuster " << decon_no_hits_gut << " nicht gefunden." <<endl;
    cout << "Decon2LS hat von " << mono_overlap.size() << " ueberlappenden Isotopenmuster " << decon_no_hits_overlap << " nicht gefunden." <<endl;

    cout << "Die Guten sind: ";
    for(int i = 0; i< decon_not_found_gut.size();++i)
    {
      cout << decon_not_found_gut[i]<<", ";
    }
    cout<< endl;

    cout << "Die ueberlappenden sind: ";
    for(int i = 0; i< decon_not_found_overlap.size();++i)
    {
      cout << decon_not_found_overlap[i]<<", ";
    }
    cout<< endl;

    cout << "OpenMS hat von " << mono_gut.size() << " guten Isotopenmuster " << open_hits_gut << " erfolgreich gefunden." <<endl;
    cout << "OpenMS hat von " << mono_overlap.size() << " ueberlappenden Isotopenmuster " << open_hits_overlap << " erfolgreich gefunden." <<endl;

    cout << "Die Guten sind: ";
    for(int i = 0; i< open_found_gut.size();++i)
    {
      cout << open_found_gut[i]<<", ";
    }
    cout<< endl;

    cout << "Die ueberlappenden sind: ";
    for(int i = 0; i< open_found_overlap.size();++i)
    {
      cout << open_found_overlap[i]<<", ";
    }
    cout<< endl;

    cout << "OpenMS hat von " << mono_gut.size() << " guten Isotopenmuster " << open_no_hits_gut << " nicht gefunden." <<endl;
    cout << "OpenMS hat von " << mono_overlap.size() << " ueberlappenden Isotopenmuster " << open_no_hits_overlap << " nicht gefunden." <<endl;

    cout << "Die Guten sind: ";
    for(int i = 0; i< open_not_found_gut.size();++i)
    {
      cout << open_not_found_gut[i]<<", ";
    }
    cout<< endl;

    cout << "Die ueberlappenden sind: ";
    for(int i = 0; i< open_not_found_overlap.size();++i)
    {
      cout << open_not_found_overlap[i]<<", ";
    }
    cout<< endl;



        ///OpenMS Deisotoper Run

        for (int i = 0; i < spectra.size(); i++)
        {

          Deisotoper::deisotopeAndSingleCharge(spectra[i],
                                               fragment_tolerance,
                                               fragment_unit_ppm,
                                               min_charge,
                                               max_charge,
                                               keep_only_deisotoped,
                                               min_isopeaks,
                                               max_isopeaks,
                                               make_single_charged,
                                               annotate_charge,
                                               annotate_iso_peak_count,
                                               use_decreasing_model,
                                               start_intensity_check,
                                               add_up_intensity);

        }

        openDeiso.setSpectra(spectra);
        openDeiso.updateRanges();

        ///Berechnung der Anzahl an monoisotopischen und nicht monoisotopischen Peaks
/*
        ///Des Centroided Datensatzes

        for(auto spec:exp_data)
        {
          all_peaks += spec.size();
          buffer_found.clear(true);
          buffer_not_found.clear(true);
          buffer_found.setRT(spec.getRT());
          buffer_not_found.setRT(spec.getRT());
          for(auto peak:spec)
          {

            bool found = false;
            auto peak_RT = spec.getRT();
            auto peak_mz = peak.getMZ();
            for(auto feature:featureMap)
            {
              double rt_min = feature.getConvexHulls()[0].getHullPoints()[0][0];
              double rt_max = feature.getConvexHulls()[0].getHullPoints()[2][0];
              double mz_min = feature.getConvexHulls()[0].getHullPoints()[0][1];
              double mz_max = feature.getConvexHulls()[0].getHullPoints()[2][1];
              if((peak_RT <= rt_max) & (peak_RT >= rt_min) & (peak_mz <= mz_max) & (peak_mz >= mz_min))
              {
                //cout << "FOUND = true Peak: " << peak << endl;
                buffer_found.push_back(peak);
                found=true;
                break;
              }
            }

            if(found == false)
            {
              //cout << "FOUND == FALSE PEAK: " << peak << endl;
              buffer_not_found.push_back(peak);
            }
          }

          mono_peaks.addSpectrum(buffer_found);
          non_mono_peaks.addSpectrum(buffer_not_found);

        }
        non_mono_peaks.updateRanges();
        mono_peaks.updateRanges();
        cout << "all " <<all_peaks << endl;
        cout << "all non " <<non_mono_peaks.getSize() << endl;
        cout << "all mono" <<mono_peaks.getSize() << endl;

          ///Des Raw Datensatzes (anteilig nur für ein bestimmtes Fenster)

          for(auto spec:exp_data_raw)
          {
            all_peaks += spec.size();
            buffer_found.clear(true);
            buffer_not_found.clear(true);
            buffer_found.setRT(spec.getRT());
            buffer_not_found.setRT(spec.getRT());
            if (spec.getRT() >= 249.177 & spec.getRT() <= 609.177)
            {
              for (auto peak:spec)
              {

                bool found = false;
                auto peak_RT = spec.getRT();
                auto peak_mz = peak.getMZ();
                for (auto feature:featureMap)
                {
                  double rt_min = feature.getConvexHulls()[0].getHullPoints()[0][0];
                  double rt_max = feature.getConvexHulls()[0].getHullPoints()[2][0];
                  double mz_min = feature.getConvexHulls()[0].getHullPoints()[0][1];
                  double mz_max = feature.getConvexHulls()[0].getHullPoints()[2][1];
                  if ((peak_RT <= rt_max) & (peak_RT >= rt_min) & (peak_mz <= mz_max) & (peak_mz >= mz_min))
                  {
                    //cout << "FOUND = true Peak: " << peak << endl;
                    buffer_found.push_back(peak);
                    found = true;
                    break;
                  }
                }

                if (found == false)
                {
                  //cout << "FOUND == FALSE PEAK: " << peak << endl;
                  buffer_not_found.push_back(peak);
                }
              }

              mono_peaks_raw.addSpectrum(buffer_found);
              non_mono_peaks_raw.addSpectrum(buffer_not_found);

            }
          }
          non_mono_peaks_raw.updateRanges();
          mono_peaks_raw.updateRanges();
          cout << "all raw " <<all_peaks_raw << endl;
          cout << "all non raw " <<non_mono_peaks_raw.getSize() << endl;
          cout << "all mono raw " <<mono_peaks_raw.getSize() << endl;
*/

        ///Berechnung der True Positiv und false Positiv Werte

/*
          ///Durchgehen des MSExperiments, dass den Hardklor Output darstellt
          for(auto spec:hardklor)
          {
            rt = spec.getRT();
            buffer_found.clear(true);
            buffer_not_found.clear(true);
            buffer_found.setRT(spec.getRT());
            buffer_not_found.setRT(spec.getRT());
            //cout << "Untersuchen des Spectrums mit RT = "<<rt << endl;
            //Durchgehen aller MSSpectren im Hardklor output
            for(auto peak:spec)
            {
              peak_found = false;
              mz = peak.getMZ();
              //cout << "Untersuchen des Peaks mit MZ= " << mz << endl;
              //Für jeden Peak in den Spectren checken ob er in der convexen Hülle eines der Features der FeatureXML steckt
              for(auto feature:featureMap)
              {

                double rt_min = feature.getConvexHulls()[0].getHullPoints()[0][0];
                double rt_max = feature.getConvexHulls()[0].getHullPoints()[2][0];
                double mz_min = feature.getConvexHulls()[0].getHullPoints()[0][1];
                double mz_max = feature.getConvexHulls()[0].getHullPoints()[2][1];

                //cout << "Untersuchen des Features mit RT " << rt_min << " - " <<rt_max << " MZ " << mz_min << " - " << mz_max << endl;

                if((rt <= rt_max) & (rt >= rt_min) & (mz <= mz_max) & (mz >= mz_min))
                {
                  buffer_found.push_back(peak);
                  peak_found = true;
                  break;
                }
                else
                {
                  if((rt <= rt_max+0.1) & (rt >= rt_min-0.1) & (mz <= mz_max+0.01) & (mz >= mz_min-0.01))
                  {
                    buffer_found.push_back(peak);
                    peak_found = true;
                    break;
                  }
                }
              }

              if(!peak_found)
              {
                buffer_not_found.push_back(peak);
              }
            }
            hardklor_true_pos.addSpectrum(buffer_found);
            hardklor_false_pos.addSpectrum(buffer_not_found);
          }

          hardklor_true_pos.updateRanges();
          hardklor_false_pos.updateRanges();

          cout << "Hardklor hat " << hardklor_true_pos.getSize() << " Peaks korrekt gefunden und " << hardklor_false_pos.getSize() << " Peaks falsch gefunden." << endl;

          ///FALSE NEGATIV ERSTELLUNG Hardklor
          for(auto mono_spec:mono_peaks)
          {
            double mono_rt = mono_spec.getRT();
            buffer_found.clear(true);
            buffer_found.setRT(mono_spec.getRT());

            for(auto mono_peak:mono_spec)
            {
              bool found = false;

              double mono_mz = mono_peak.getMZ();
              double mono_int = mono_peak.getIntensity();
              //cout << "mono peak: rt " <<mono_rt << " mz " << mono_mz <<endl;

              for (auto spec:hardklor_true_pos)
              {
                rt =spec.getRT();

                if(rt == mono_rt)
                {
                  for(int i = 0; i<spec.size();++i)
                  {
                    double peak_mz = spec[i].getMZ();
                    //cout << "Peak mz " << peak_mz << endl;
                    //cout << "ABS "<<abs(peak_mz-mono_mz)<<endl;
                    if( abs(mono_mz - peak_mz) < 0.0001 )
                    {
                      //cout << "GEFUNDEN!"<< endl;
                      found = true;
                      break;
                    }
                  }
                }
                if(found)
                {
                  break;
                }
              }
              if(!found)
              {
                buffer_found.push_back(mono_peak);
              }
            }
            hardklor_false_neg.addSpectrum(buffer_found);
          }

          hardklor_false_neg.updateRanges();

          ofstream output;
          output.open("/buffer/ag_bsc/pmsb_2021/Rieckert/data/hardklor_false_negativ.txt");

          for(auto spec:hardklor_false_neg)
          {
            rt = spec.getRT();
            for(auto peak:spec)
            {
              output << "Peak: RT: " << rt << " MZ: " << peak.getMZ() << " Int: " << peak.getIntensity()<< "\n" << endl;
            }
          }
          output.close();

          cout << hardklor_false_neg.getSize() << endl;

        ///Durchgehen des MSExperiments, dass den Decon2LS Output darstellt

        for(auto spec:decon)
        {
          rt = spec.getRT();
          buffer_found.clear(true);
          buffer_not_found.clear(true);
          buffer_found.setRT(spec.getRT());
          buffer_not_found.setRT(spec.getRT());
          //Durchgehen aller MSSpectren im Decon2LS output
          for(auto peak:spec)
          {
            peak_found = false;
            mz = peak.getMZ();
            //Für jeden Peak in den Spectren checken ob er in der convexen Hülle eines der Features der FeatureXML steckt
            for(auto feature:featureMap)
            {

              double rt_min = feature.getConvexHulls()[0].getHullPoints()[0][0];
              double rt_max = feature.getConvexHulls()[0].getHullPoints()[2][0];
              double mz_min = feature.getConvexHulls()[0].getHullPoints()[0][1];
              double mz_max = feature.getConvexHulls()[0].getHullPoints()[2][1];

              if((rt <= rt_max) & (rt >= rt_min) & (mz <= mz_max) & (mz >= mz_min))
              {
                buffer_found.push_back(peak);
                peak_found = true;
                break;
              }
              else
              {
                if((rt <= rt_max+0.1) & (rt >= rt_min-0.1) & (mz <= mz_max+0.01) & (mz >= mz_min-0.01))
                {
                  buffer_found.push_back(peak);
                  peak_found = true;
                  break;
                }
              }
            }

            if(!peak_found)
            {
              buffer_not_found.push_back(peak);
            }
          }
          decon_true_pos.addSpectrum(buffer_found);
          decon_false_pos.addSpectrum(buffer_not_found);
        }

        decon_true_pos.updateRanges();
        decon_false_pos.updateRanges();
        cout << "Decon2LS hat " << decon_true_pos.getSize() << " Peaks korrekt gefunden und " << decon_false_pos.getSize() << " Peaks falsch gefunden." << endl;


        //Durchgehen des MSExperiments, dass den OpenMS Deisotoper Output darstellt
        for(auto spec:openDeiso)
        {
          rt = spec.getRT();
          buffer_found.clear(true);
          buffer_not_found.clear(true);
          buffer_found.setRT(spec.getRT());
          buffer_not_found.setRT(spec.getRT());
          //Durchgehen aller MSSpectren im Hardklor output
          for(auto peak:spec)
          {
            peak_found = false;
            mz = peak.getMZ();
            //Für jeden Peak in den Spectren checken ob er in der convexen Hülle eines der Features der FeatureXML steckt
            for(auto feature:featureMap)
            {

              double rt_min = feature.getConvexHulls()[0].getHullPoints()[0][0];
              double rt_max = feature.getConvexHulls()[0].getHullPoints()[2][0];
              double mz_min = feature.getConvexHulls()[0].getHullPoints()[0][1];
              double mz_max = feature.getConvexHulls()[0].getHullPoints()[2][1];

              if((rt <= rt_max) & (rt >= rt_min) & (mz <= mz_max) & (mz >= mz_min))
              {
                buffer_found.push_back(peak);
                peak_found = true;
                break;
              }
              else
              {
                if((rt <= rt_max+0.1) & (rt >= rt_min-0.1) & (mz <= mz_max+0.01) & (mz >= mz_min-0.01))
                {
                  buffer_found.push_back(peak);
                  peak_found = true;
                  break;
                }
              }
            }

            if(!peak_found)
            {
              buffer_not_found.push_back(peak);
            }
          }
          openDeiso_true_pos.addSpectrum(buffer_found);
          openDeiso_false_pos.addSpectrum(buffer_not_found);
        }

        openDeiso_true_pos.updateRanges();
        openDeiso_false_pos.updateRanges();

        cout << "Der OpenMS Deisotoper hat " << openDeiso_true_pos.getSize() << " Peaks korrekt gefunden und " << openDeiso_false_pos.getSize() << " Peaks falsch gefunden." << endl;

        */

    return;
  }


  ///Veraltete Version der Kullback-Leibler-Divergenz berechnung, ist in realKLD kompakter Implementiert
  void DeisotoperRieckert::kldVerify()
  {
     String input = "/buffer/ag_bsc/pmsb_2021/Rieckert/data/KLD-protein-samples.FASTA";

     ProteaseDigestion digestor;

     //FineIsotopePatternGenerator generator;
     CoarseIsotopePatternGenerator generator;

     FASTAFile file;
     vector<FASTAFile::FASTAEntry> proteins;

     vector<AASequence> undigested_seq;
     vector<vector<AASequence>> digested_seq;
     vector<AASequence> tmpAASVec;
     AASequence tmpAAS;

     vector<EmpiricalFormula> tmpEmpiVec;
     vector<vector<EmpiricalFormula>> digested_formulas;
     EmpiricalFormula tmp_formula;

     vector<IsotopeDistribution> tmp_isotopes;
     vector<vector<IsotopeDistribution>> digested_distributions;
     IsotopeDistribution tmp_distribution;

     IsotopeModel averagine_generator;

     vector<EmpiricalFormula> tmp_averagine_formula_vec;
     vector<vector<EmpiricalFormula>> averagine_formulas;
     EmpiricalFormula tmp_averagine_formula;

     vector<IsotopeDistribution> tmp_averagine_distribution_vec;
     vector<vector<IsotopeDistribution>> averagine_distributions;
     IsotopeDistribution tmp_averagine_distribution;

     vector<double> tmp_kld_vec;
     vector<vector<double>> kld_values;
     double tmp_kld;

     vector<double> kld_n2;
     vector<double> kld_n3;
     vector<double> kld_n4;
     vector<double> kld_n5;
     vector<double> kld_n6;
     vector<double> kld_n7;
     vector<double> kld_n8;
     vector<double> kld_n9;


     ///read FASTA

     file.load(input,proteins);

     ///convert FASTAEntrys into AASequenzes

     for(auto fasta:proteins)
     {
       tmpAAS = AASequence::fromString(fasta.sequence);
       undigested_seq.push_back(tmpAAS);
     }

     cout << "Undigested size: " << undigested_seq.size() << endl;

     ///Digest the proteins using Trypsin

     digestor.setEnzyme("trypsin");
     digestor.setMissedCleavages(2);

     for(auto prot:undigested_seq)
     {
       digestor.digest(prot,tmpAASVec);
       digested_seq.push_back(tmpAASVec);
       tmpAASVec.clear();
     }

     cout << "Digested size: " << digested_seq.size() << endl;

     /// generating the empirical formula for the peptides which are used to generate the theoretical isotope distributions

     for(auto pep_vec:digested_seq)
     {
       for(auto pep:pep_vec)
       {
         tmp_formula = pep.getFormula(Residue::Full);
         tmpEmpiVec.push_back(tmp_formula);
       }
       digested_formulas.push_back(tmpEmpiVec);
       tmpEmpiVec.clear();
     }

     for(auto form_vec:digested_formulas)
     {
       for(auto formula:form_vec)
       {
         tmp_distribution = generator.run(formula);
         tmp_isotopes.push_back(tmp_distribution);
       }
       digested_distributions.push_back(tmp_isotopes);
       tmp_isotopes.clear();
     }

      //cout << digested_seq[0][0].getAverageWeight(Residue::Internal,1) << " : " << digested_distributions[0][0].averageMass() << endl;
     ///Generating a vector filled with the empirical formulas of the averagine models using the average mass of the digested peptide distributions

     for(auto pep_vec:digested_seq) //digested_distributions
     {
       for(auto pep:pep_vec)
       {
         Param mean ;
         mean.setValue("statistics:mean", pep.getAverageWeight(Residue::Full,1), "Centroid m/z (as opposed to monoisotopic m/z).", ListUtils::create<String>("advanced"));
         averagine_generator.setParameters(mean);
         tmp_averagine_formula = averagine_generator.getFormula();
         tmp_averagine_formula_vec.push_back(tmp_averagine_formula);
       }
       averagine_formulas.push_back(tmp_averagine_formula_vec);
       tmp_averagine_formula_vec.clear();
     }

     ///Using these averagine formulas to generate the averagine isotopic distributions

     for(auto form_vec:averagine_formulas)
     {
       for(auto formula:form_vec)
       {
         averagine_generator.setSamples(formula);
         tmp_averagine_distribution = averagine_generator.getIsotopeDistribution();
         tmp_averagine_distribution_vec.push_back(tmp_averagine_distribution);
       }
       averagine_distributions.push_back(tmp_averagine_distribution_vec);
       tmp_averagine_distribution_vec.clear();
     }

     /// Application of the Kullback-Leibler-divergence

     // KL(P,A) = sum(P(x) * log(P(x)/A(x))) over all x || P = peptide distribution; A = Averagine distribution ; x = m/z of the n-th Peak of the distribution

     ///N=2

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=2 & tmp_dist_ave.size()>=2)
         {
           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity();
           double tmp_ave = tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity();
           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           //cout << "tmp_dist d1 " << tmp_d_1 << " d2 " << tmp_d_2 << " a1 " << tmp_a_1 << " a2 " << tmp_a_2 << endl;
           kld_n2.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(2, tmp_dist, tmp_dist_ave));

         }

        // kld_n2.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(2,digested_distributions[d][i],tmp_averagine_distribution_vec[i]));
       }
     }
     kld_n2.erase(remove(kld_n2.begin(),kld_n2.end(), -9999),kld_n2.end());

     ///N=3

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=3 & tmp_dist_ave.size()>=3)
         {

           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);

           kld_n3.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(3, tmp_dist, tmp_dist_ave));
         }
       }
     }
     kld_n3.erase(remove(kld_n3.begin(),kld_n3.end(), -9999),kld_n3.end());

     ///N=4

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=4 & tmp_dist_ave.size()>=4)
         {
           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                               tmp_dist[3].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
               tmp_dist_ave[3].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
           double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
           double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist[3].setIntensity(tmp_d_4);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);
           tmp_dist_ave[3].setIntensity(tmp_a_4);

           kld_n4.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(4, tmp_dist, tmp_dist_ave));
         }
       }
     }
     kld_n4.erase(remove(kld_n4.begin(),kld_n4.end(), -9999),kld_n4.end());

     ///N=5

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=5 & tmp_dist_ave.size()>=5)
         {

           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                               tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
               tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
           double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
           double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
           double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
           double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist[3].setIntensity(tmp_d_4);
           tmp_dist[4].setIntensity(tmp_d_5);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);
           tmp_dist_ave[3].setIntensity(tmp_a_4);
           tmp_dist_ave[4].setIntensity(tmp_a_5);

           kld_n5.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(5, tmp_dist, tmp_dist_ave));
         }
       }
     }

     kld_n5.erase(remove(kld_n5.begin(),kld_n5.end(), -9999),kld_n5.end());

     ///N=6

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=6 & tmp_dist_ave.size()>=6)
         {

           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                               tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
               tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
           double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
           double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
           double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
           double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
           double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
           double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist[3].setIntensity(tmp_d_4);
           tmp_dist[4].setIntensity(tmp_d_5);
           tmp_dist[5].setIntensity(tmp_d_6);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);
           tmp_dist_ave[3].setIntensity(tmp_a_4);
           tmp_dist_ave[4].setIntensity(tmp_a_5);
           tmp_dist_ave[5].setIntensity(tmp_a_6);

           kld_n6.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(6, tmp_dist, tmp_dist_ave));
         }
       }
     }

     kld_n6.erase(remove(kld_n6.begin(),kld_n6.end(), -9999),kld_n6.end());

     ///N=7

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=7 & tmp_dist_ave.size()>=7)
         {

           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                               tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity() + tmp_dist[6].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
               tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity() + tmp_dist_ave[6].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
           double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
           double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
           double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;
           double tmp_d_7 = tmp_dist[6].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
           double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
           double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
           double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;
           double tmp_a_7 = tmp_dist_ave[6].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist[3].setIntensity(tmp_d_4);
           tmp_dist[4].setIntensity(tmp_d_5);
           tmp_dist[5].setIntensity(tmp_d_6);
           tmp_dist[6].setIntensity(tmp_d_7);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);
           tmp_dist_ave[3].setIntensity(tmp_a_4);
           tmp_dist_ave[4].setIntensity(tmp_a_5);
           tmp_dist_ave[5].setIntensity(tmp_a_6);
           tmp_dist_ave[6].setIntensity(tmp_a_7);

           kld_n7.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(7, tmp_dist, tmp_dist_ave));
         }
       }
     }

     kld_n7.erase(remove(kld_n7.begin(),kld_n7.end(), -9999),kld_n7.end());

     ///N=8

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=8 & tmp_dist_ave.size()>=8)
         {

           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                               tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity() + tmp_dist[6].getIntensity() + tmp_dist[7].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
               tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity() + tmp_dist_ave[6].getIntensity() + tmp_dist_ave[7].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
           double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
           double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
           double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;
           double tmp_d_7 = tmp_dist[6].getIntensity() / tmp_digest;
           double tmp_d_8 = tmp_dist[7].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
           double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
           double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
           double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;
           double tmp_a_7 = tmp_dist_ave[6].getIntensity() / tmp_ave;
           double tmp_a_8 = tmp_dist_ave[7].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist[3].setIntensity(tmp_d_4);
           tmp_dist[4].setIntensity(tmp_d_5);
           tmp_dist[5].setIntensity(tmp_d_6);
           tmp_dist[6].setIntensity(tmp_d_7);
           tmp_dist[7].setIntensity(tmp_d_8);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);
           tmp_dist_ave[3].setIntensity(tmp_a_4);
           tmp_dist_ave[4].setIntensity(tmp_a_5);
           tmp_dist_ave[5].setIntensity(tmp_a_6);
           tmp_dist_ave[6].setIntensity(tmp_a_7);
           tmp_dist_ave[7].setIntensity(tmp_a_8);

           kld_n8.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(8, tmp_dist, tmp_dist_ave));
         }


       }
     }

     kld_n8.erase(remove(kld_n8.begin(),kld_n8.end(), -9999),kld_n8.end());

     ///N=9

     for(int d = 0; d< digested_distributions.size();++d)
     {
       tmp_averagine_distribution_vec.clear();
       tmp_averagine_distribution_vec = averagine_distributions[d];
       for(int i = 0; i < digested_distributions[d].size(); ++i)
       {
         digested_distributions[d][i].trimIntensities(0.001);

         IsotopeDistribution tmp_dist = digested_distributions[d][i];
         IsotopeDistribution tmp_dist_ave = tmp_averagine_distribution_vec[i];

         if(tmp_dist.size()>=9 & tmp_dist_ave.size()>=9)
         {

           double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                               tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity() + tmp_dist[6].getIntensity() + tmp_dist[7].getIntensity() + tmp_dist[8].getIntensity();
           double tmp_ave =
               tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
               tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity() + tmp_dist_ave[6].getIntensity() + tmp_dist_ave[7].getIntensity()  + tmp_dist_ave[8].getIntensity();


           double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
           double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
           double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
           double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
           double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
           double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;
           double tmp_d_7 = tmp_dist[6].getIntensity() / tmp_digest;
           double tmp_d_8 = tmp_dist[7].getIntensity() / tmp_digest;
           double tmp_d_9 = tmp_dist[8].getIntensity() / tmp_digest;

           double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
           double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
           double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
           double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
           double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
           double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;
           double tmp_a_7 = tmp_dist_ave[6].getIntensity() / tmp_ave;
           double tmp_a_8 = tmp_dist_ave[7].getIntensity() / tmp_ave;
           double tmp_a_9 = tmp_dist_ave[8].getIntensity() / tmp_ave;

           tmp_dist[0].setIntensity(tmp_d_1);
           tmp_dist[1].setIntensity(tmp_d_2);
           tmp_dist[2].setIntensity(tmp_d_3);
           tmp_dist[3].setIntensity(tmp_d_4);
           tmp_dist[4].setIntensity(tmp_d_5);
           tmp_dist[5].setIntensity(tmp_d_6);
           tmp_dist[6].setIntensity(tmp_d_7);
           tmp_dist[7].setIntensity(tmp_d_8);
           tmp_dist[8].setIntensity(tmp_d_9);
           tmp_dist_ave[0].setIntensity(tmp_a_1);
           tmp_dist_ave[1].setIntensity(tmp_a_2);
           tmp_dist_ave[2].setIntensity(tmp_a_3);
           tmp_dist_ave[3].setIntensity(tmp_a_4);
           tmp_dist_ave[4].setIntensity(tmp_a_5);
           tmp_dist_ave[5].setIntensity(tmp_a_6);
           tmp_dist_ave[6].setIntensity(tmp_a_7);
           tmp_dist_ave[7].setIntensity(tmp_a_8);
           tmp_dist_ave[8].setIntensity(tmp_a_9);


           kld_n9.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(9, tmp_dist, tmp_dist_ave));
         }


       }
     }

     kld_n9.erase(remove(kld_n9.begin(),kld_n9.end(), -9999),kld_n9.end());

     cout << "N2:" <<kld_n2.size() <<endl;

    auto itA2 = max_element(std::begin(kld_n2), std::end(kld_n2));
    cout << *itA2 <<endl;

    auto itB2 = min_element(std::begin(kld_n2), std::end(kld_n2));
    cout << *itB2 <<endl;

    cout << "N3:" <<kld_n3.size() <<endl;

    auto itA3 = max_element(std::begin(kld_n3), std::end(kld_n3));
    cout << *itA3 <<endl;

    auto itB3 = min_element(std::begin(kld_n3), std::end(kld_n3));
    cout << *itB3 <<endl;

    cout << "N4:" <<kld_n4.size() <<endl;

    auto itA4 = max_element(std::begin(kld_n4), std::end(kld_n4));
    cout << *itA4 <<endl;

    auto itB4 = min_element(std::begin(kld_n4), std::end(kld_n4));
    cout << *itB4 <<endl;

    cout << "N5:" <<kld_n5.size() <<endl;

    auto itA5 = max_element(std::begin(kld_n5), std::end(kld_n5));
    cout << *itA5 <<endl;

    auto itB5 = min_element(std::begin(kld_n5), std::end(kld_n5));
    cout << *itB5 <<endl;

    cout << "N6:" <<kld_n6.size() <<endl;

    auto itA6 = max_element(std::begin(kld_n6), std::end(kld_n6));
    cout << *itA6 <<endl;

    auto itB6 = min_element(std::begin(kld_n6), std::end(kld_n6));
    cout << *itB6 <<endl;

    cout << "N7:" <<kld_n7.size() <<endl;

    auto itA7 = max_element(std::begin(kld_n7), std::end(kld_n7));
    cout << *itA7 <<endl;

    auto itB7 = min_element(std::begin(kld_n7), std::end(kld_n7));
    cout << *itB7 <<endl;

    cout << "N8:" <<kld_n8.size() <<endl;

    auto itA8 = max_element(std::begin(kld_n8), std::end(kld_n8));
    cout << *itA8 <<endl;

    auto itB8 = min_element(std::begin(kld_n8), std::end(kld_n8));
    cout << *itB8 <<endl;

    cout << "N9:" <<kld_n9.size() <<endl;

    auto itA9 = max_element(std::begin(kld_n9), std::end(kld_n9));
    cout << *itA9 <<endl;

    auto itB9 = min_element(std::begin(kld_n9), std::end(kld_n9));
    cout << *itB9 <<endl;



    //cout << 0.0012183 * log(0.0012183/0.0012183) <<endl;




     //sort(kld_n2.begin(),kld_n2.end());
/*
                ofstream n2stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n2_normalized_mc2.txt");
                copy(kld_n2.begin(), kld_n2.end(), std::ostream_iterator<double>(n2stream, "\n"));

                ofstream n3stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n3_normalized_mc2.txt");
                copy(kld_n3.begin(), kld_n3.end(), std::ostream_iterator<double>(n3stream, "\n"));

                ofstream n4stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n4_normalized_mc2.txt");
                copy(kld_n4.begin(), kld_n4.end(), std::ostream_iterator<double>(n4stream, "\n"));

                ofstream n5stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n5_normalized_mc2.txt");
                copy(kld_n5.begin(), kld_n5.end(), std::ostream_iterator<double>(n5stream, "\n"));

                ofstream n6stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n6_normalized_mc2.txt");
                copy(kld_n6.begin(), kld_n6.end(), std::ostream_iterator<double>(n6stream, "\n"));

                ofstream n7stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n7_normalized_mc2.txt");
                copy(kld_n7.begin(), kld_n7.end(), std::ostream_iterator<double>(n7stream, "\n"));

                ofstream n8stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n8_normalized_mc2.txt");
                copy(kld_n8.begin(), kld_n8.end(), std::ostream_iterator<double>(n8stream, "\n"));

                ofstream n9stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n9_normalized_mc2.txt");
                copy(kld_n9.begin(), kld_n9.end(), std::ostream_iterator<double>(n9stream, "\n"));
*/

    /*
                     ///Testprints

                     cout << "N=2 : ";
                     for(int i = 0; i< kld_n2.size(); ++i)
                     {
                       cout << kld_n2[i] << " , ";
                     }
                     cout << endl;

                    cout << "N=3 : ";
                    for(int i = 0; i< kld_n3.size(); ++i)
                    {
                      cout << kld_n3[i] << " , ";
                    }
                    cout << endl;

                    cout << "N=4 : ";
                    for(int i = 0; i< kld_n4.size(); ++i)
                    {
                      cout << kld_n4[i] << " , ";
                    }
                    cout << endl;

                    cout << "N=5 : ";
                    for(int i = 0; i< kld_n5.size(); ++i)
                    {
                      cout << kld_n5[i] << " , ";
                    }
                    cout << endl;

                    cout << "N=6 : ";
                    for(int i = 0; i< kld_n6.size(); ++i)
                    {
                      cout << kld_n6[i] << " , ";
                    }
                    cout << endl;

                    cout << "N=7 : ";
                    for(int i = 0; i< kld_n7.size(); ++i)
                    {
                      cout << kld_n7[i] << " , ";
                    }
                    cout << endl;

                    cout << "N=8 : ";
                    for(int i = 0; i< kld_n8.size(); ++i)
                    {
                      cout << kld_n8[i] << " , ";
                    }
                    cout << endl;

                    cout << "N=9 : ";
                    for(int i = 0; i< kld_n9.size(); ++i)
                    {
                      cout << kld_n9[i] << " , ";
                    }
                    cout << endl;

                       for(auto dist_vec:digested_distributions)
    {
      for(auto dist:dist_vec)
      {
        for(auto peak:dist)
        {
          cout << peak << endl;
        }
        break;
      }
      break;
    }

    for(auto dist_vec:averagine_distributions)
    {
      for(auto dist:dist_vec)
      {
        for(auto peak:dist)
        {
          cout << peak << endl;
        }
        break;
      }
      break;
    }

                */

  }

  ///Berechnet die Kullback-Leibler-Divergenz für zwei Isotopendistributionen
  double DeisotoperRieckert::kullbackLeiblerDivergence(int n, IsotopeDistribution pep_distri, IsotopeDistribution ave_distri)
  {
    if(n > pep_distri.size() | n > ave_distri.size())
    {
      ///Fehlerwert anhand dessen ich die Ausgabe filterte
      return -9999;
    }

    double kld = 0;

    for(int i = 0; i< n; ++i)
    {
      kld += pep_distri[i].getIntensity() * log((pep_distri[i].getIntensity() / ave_distri[i].getIntensity()));
    }
    return kld;
  }

  ///Überprüft ob die Peaks einer Isotopendistribution die korrekten Intensitätsverhältnisse haben
  bool DeisotoperRieckert::distributionValidation(OpenMS::IsotopeDistribution input)
  {
    double start = input[0].getIntensity();
    bool is_valid = false;


    ///Fall 1: Alle sind kleiner als der erste.
    if(input[1].getIntensity() < start)
    {
      ///Wenn es mehr als 2 Peaks in der Distribution gibt überprüfe die restlichen, falls nicht return true
      if(input.size()>2)
      {
        for (int i = 2; i < input.size(); ++i)
        {
          /// Check ob der Peak kleiner ist als sein Vorgänger
          if(input[i].getIntensity() < input[i-1].getIntensity())
          {
            ///  Falls wir im letzen peak angekommen sind dann return true;
            if(i==input.size()-1)
            {
              return true;
            }
            continue;
          }
          ///Falls der aktuelle Peak größer als der Vorgänger ist return false;
          else
          {
            return false;
          }
        }
      }
      else
      {
        return true;
      }
    }

    ///Fall 2: Alle sind größer als der erste (evt gibt es solche Distriibutionen gar nicht)
    ///Fall 3: Die Intensity Werte steigen erst und fallen dann wieder (Glockenform)

    if(input[1].getIntensity() > start)
    {
      ///Wenn es mehr als 2 Peaks in der Distribution gibt überprüfe die restlichen, falls nicht return true
      if(input.size()>2)
      {
        for (int i = 2; i < input.size(); ++i)
        {
          /// Check ob der Peak größer ist als sein Vorgänger
          if(input[i].getIntensity() > input[i-1].getIntensity())
          {
            ///  Falls wir im letzen peak angekommen sind dann return true;
            if(i==input.size()-1)
            {
              return true;
            }
            continue;
          }
            ///Falls der aktuelle Peak kleiner als der Vorgänger ist überprüfe ob die nachfolgenden absteigend sind
          else
          {
            ///Falls wir noch nicht im letzten angekommen waren
            if(i!=input.size()-1)
            {
              for (int j = i+1; j < input.size(); ++j)
              {
                /// Check ob der Peak kleiner ist als sein Vorgänger
                if(input[j].getIntensity() < input[j-1].getIntensity())
                {
                  ///Falls wir im letzten angekommen sind return true
                  if(j==input.size())
                  {
                    return true;
                  }
                  continue;
                }
                ///Fall der Peak größer ist als sein Vorgänger haben wir den Fall dass die werte erst aufsteigen dann fallen und dann wieder aufsteigen!! -> return false
                else
                {
                  return false;
                }
              }
            }
            ///Falls wir beim letzten waren und dieser als einziger kleiner war return true
            else
            {
              return true;
            }
          }
        }
      }
      else
      {
        return true;
      }
    }

    return is_valid;
  }

  ///Erstellt für eine gegebene Masse die entsprechende Averagineisotopendistribution
  OpenMS::IsotopeDistribution DeisotoperRieckert::averagineDistribuitonGenerator(double mass)
  {
    IsotopeModel averagine_generator;
    EmpiricalFormula averagine_formula;

    Param mean;
    mean.setValue("statistics:mean", mass, "Centroid m/z (as opposed to monoisotopic m/z).", ListUtils::create<String>("advanced"));
    averagine_generator.setParameters(mean);
    averagine_formula = averagine_generator.getFormula();
    averagine_generator.setSamples(averagine_formula);

    return averagine_generator.getIsotopeDistribution();

  }

  ///Berechnung der Kullback-Leibler-Divergenzen für den Testdatensatz
  void DeisotoperRieckert::realKLD()
  {
    String inputMZ = "/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/trimmed/MC_bud3_Trypsin_1_trimmed_1412_92.mzML";

    MzMLFile mzFile;

    MSExperiment exp_raw;
    MSExperiment exp_picked;
    MSExperiment exp_deiso;

    double fragment_tolerance = 5.0;
    bool fragment_unit_ppm = true;
    int min_charge = 2;
    int max_charge = 4;
    bool keep_only_deisotoped = true;
    unsigned int min_isopeaks = 2;
    unsigned int max_isopeaks = 5;
    bool make_single_charged = false;
    bool annotate_charge = true;
    bool annotate_iso_peak_count = true;
    bool use_decreasing_model = true;
    unsigned int start_intensity_check = 1;
    bool add_up_intensity = false;

    vector<double> kld_n2;
    vector<double> kld_n3;
    vector<double> kld_n4;
    vector<double> kld_n5;
    vector<double> kld_n6;
    vector<double> kld_n7;
    vector<double> kld_n8;
    vector<double> kld_n9;

    vector<IsotopeDistribution> n2_dist_tresh_hits;
    vector<IsotopeDistribution> n2_ave_tresh_hits;
    vector<IsotopeDistribution> n3_dist_tresh_hits;
    vector<IsotopeDistribution> n3_ave_tresh_hits;
    vector<IsotopeDistribution> n4_dist_tresh_hits;
    vector<IsotopeDistribution> n4_ave_tresh_hits;
    vector<IsotopeDistribution> n5_dist_tresh_hits;
    vector<IsotopeDistribution> n5_ave_tresh_hits;
    vector<IsotopeDistribution> n6_dist_tresh_hits;
    vector<IsotopeDistribution> n6_ave_tresh_hits;

    vector<int> chargeStates = {2,3,4};
    vector<Peak1D> mono_peaks;



    mzFile.load(inputMZ,exp_raw);
    PeakPickerHiRes pp;
    pp.pickExperiment(exp_raw,exp_picked);


    vector<MSSpectrum> spectra = exp_picked.getSpectra();

    for(auto spec:spectra)
    {
      Deisotoper::deisotopeAndSingleCharge(spec,
                                           fragment_tolerance,
                                           fragment_unit_ppm,
                                           min_charge,
                                           max_charge,
                                           keep_only_deisotoped,
                                           min_isopeaks,
                                           max_isopeaks,
                                           make_single_charged,
                                           annotate_charge,
                                           annotate_iso_peak_count,
                                           use_decreasing_model,
                                           start_intensity_check,
                                           add_up_intensity);
    }

    exp_deiso.setSpectra(spectra);
    exp_deiso.updateRanges();

    for(auto spec:exp_deiso)
    {
      for(auto peak:spec)
      {
        mono_peaks.push_back(peak);
      }
    }

    /// Erstellen der Distributionen der Monoisotopischen Peaks die OpenMS gefunden hat

    vector<IsotopeDistribution> isotope_distributions;
    vector<IsotopeDistribution> averagine_distributions;

    for(auto mono_peak:mono_peaks)
    {
      IsotopeDistribution tmp_distri;
      vector<IsotopeDistribution> potential_distributions;
      IsotopeDistribution best_candidate;
      double mono_mz = mono_peak.getMZ();
      vector<int> charge_hits;

      for(auto charge:chargeStates)
      {
        vector<Peak1D> charged_distri;
        charged_distri.push_back(mono_peak);


        for(auto ori_spec:exp_picked)
        {
          int count = 1;
          for(auto ori_peak:ori_spec)
          {
            double ori_mz = ori_peak.getMZ();
            if( (abs( ori_mz - (mono_mz + ((1.0/charge)*count)) ) <= (0.015)) )
            {
              charged_distri.push_back(ori_peak);
              count +=1;
            }
            else
            {
              if(( (abs( ori_mz - (mono_mz - ((1.0/charge)) ) ) <= 0.01) & (ori_peak.getIntensity() < mono_peak.getIntensity()) ) )
              {
                charged_distri.push_back(ori_peak);
              }
            }
          }
        }

        if(charged_distri.size() >1 ) //Wenn mehr als nur der monoisotopische Peak in der potentiellen Distribution enthalten ist.
        {
          //IsotopeDistribution tmp_distri;
          tmp_distri.set(charged_distri);
          tmp_distri.sortByMass();

          if(DeisotoperRieckert::distributionValidation(tmp_distri))
          {
            potential_distributions.push_back(tmp_distri);
            charge_hits.push_back(charge);
          }
        }

        charged_distri.clear();

      }

      if(potential_distributions.size() > 1)
      {
        //cout << "MEHR ALS EINE POTENTIELLE!!!!! " << mono_mz<<endl;
        //cout << potential_distributions.size() <<endl;
        //cout << charge_hits.size() << endl;

        for(int i = 0 ; i <potential_distributions.size(); ++i)
        {
          ///Da die Anzahl an genau passenden Isotopendistributionen für meine Arbeit ausreichte ignorierte ich alle Isotopendistributionen, bei denen es zu mehr als einem Treffer für die Abstände zum monoisotopischen Peak kam.
          continue;
          /* Hier muss nun die entscheidung getroffen werden Welche der potentiellen Distributionen die bessere ist
           * Dann muss potential_distributions[i] in isotope_distributions gepusht werden
           * und averagineDistribuitongenerator(mono_mz*charge_hits[i]) in averagine_distributions
           */
        }
      }
      else{

        if(potential_distributions.size() == 0)
        {
          //cout << "KEINE POTENTIELLE " << mono_mz << endl;
          continue;
        }
        else
        {
          //cout << "Genau eine pot. " <<endl;
          isotope_distributions.push_back(potential_distributions[0]);
          averagine_distributions.push_back(averagineDistribuitonGenerator((mono_mz*charge_hits[0])));
        }
      }
      potential_distributions.clear();
      best_candidate.clear();
      tmp_distri.clear();
      charge_hits.clear();

    }

    cout << isotope_distributions.size() << endl;
    cout << averagine_distributions.size() << endl;


    ///KL-D Berechnung

    ///N=2
    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=2 & tmp_dist_ave.size()>=2)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity();
        double tmp_ave = tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity();
        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);

        double tmp_kld = DeisotoperRieckert::kullbackLeiblerDivergence(2, tmp_dist, tmp_dist_ave);

        kld_n2.push_back(tmp_kld);

        if(tmp_kld <= 0.05 & tmp_kld >= 0.045)
        {
          n2_ave_tresh_hits.push_back(tmp_dist_ave);
          n2_dist_tresh_hits.push_back(tmp_dist);



        }


      }
    }
    kld_n2.erase(remove(kld_n2.begin(),kld_n2.end(), -9999),kld_n2.end());

    ///N=3

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=3 & tmp_dist_ave.size()>=3)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity();
        double tmp_ave = tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);

        double tmp_kld = DeisotoperRieckert::kullbackLeiblerDivergence(3, tmp_dist, tmp_dist_ave);

        kld_n3.push_back(tmp_kld);

        if(tmp_kld <= 0.1 & tmp_kld >= 0.09)
        {
          n3_ave_tresh_hits.push_back(tmp_dist_ave);
          n3_dist_tresh_hits.push_back(tmp_dist);



        }

        //kld_n3.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(3, tmp_dist, tmp_dist_ave));

      }
    }
    kld_n3.erase(remove(kld_n3.begin(),kld_n3.end(), -9999),kld_n3.end());

    ///N=4

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=4 & tmp_dist_ave.size()>=4)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                            tmp_dist[3].getIntensity();
        double tmp_ave =
            tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
            tmp_dist_ave[3].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
        double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
        double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist[3].setIntensity(tmp_d_4);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);
        tmp_dist_ave[3].setIntensity(tmp_a_4);

        double tmp_kld = DeisotoperRieckert::kullbackLeiblerDivergence(4, tmp_dist, tmp_dist_ave);

        kld_n4.push_back(tmp_kld);

        if(tmp_kld <= 0.2 & tmp_kld >= 0.15)
        {
          n4_ave_tresh_hits.push_back(tmp_dist_ave);
          n4_dist_tresh_hits.push_back(tmp_dist);

        }

        //kld_n4.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(4, tmp_dist, tmp_dist_ave));

      }
    }
    kld_n4.erase(remove(kld_n4.begin(),kld_n4.end(), -9999),kld_n4.end());

    ///N=5

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=5 & tmp_dist_ave.size()>=5)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                            tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity();
        double tmp_ave =
            tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
            tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
        double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
        double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
        double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
        double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist[3].setIntensity(tmp_d_4);
        tmp_dist[4].setIntensity(tmp_d_5);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);
        tmp_dist_ave[3].setIntensity(tmp_a_4);
        tmp_dist_ave[4].setIntensity(tmp_a_5);

        double tmp_kld = DeisotoperRieckert::kullbackLeiblerDivergence(5, tmp_dist, tmp_dist_ave);

        kld_n5.push_back(tmp_kld);

        if(tmp_kld <= 0.4 & tmp_kld >= 0.35)
        {
          n5_ave_tresh_hits.push_back(tmp_dist_ave);
          n5_dist_tresh_hits.push_back(tmp_dist);


        }


        //kld_n5.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(5, tmp_dist, tmp_dist_ave));
      }
    }
    kld_n5.erase(remove(kld_n5.begin(),kld_n5.end(), -9999),kld_n5.end());

    ///N=6

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=6 & tmp_dist_ave.size()>=6)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                            tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity();
        double tmp_ave =
            tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
            tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
        double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
        double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
        double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
        double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
        double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
        double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist[3].setIntensity(tmp_d_4);
        tmp_dist[4].setIntensity(tmp_d_5);
        tmp_dist[5].setIntensity(tmp_d_6);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);
        tmp_dist_ave[3].setIntensity(tmp_a_4);
        tmp_dist_ave[4].setIntensity(tmp_a_5);
        tmp_dist_ave[5].setIntensity(tmp_a_6);

        double tmp_kld = DeisotoperRieckert::kullbackLeiblerDivergence(6, tmp_dist, tmp_dist_ave);

        kld_n6.push_back(tmp_kld);

        if(tmp_kld <= 0.6 & tmp_kld >= 0.55)
        {
          n6_ave_tresh_hits.push_back(tmp_dist_ave);
          n6_dist_tresh_hits.push_back(tmp_dist);


        }

       // kld_n6.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(6, tmp_dist, tmp_dist_ave));
      }
    }
    kld_n6.erase(remove(kld_n6.begin(),kld_n6.end(), -9999),kld_n6.end());

    ///N=7

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=7 & tmp_dist_ave.size()>=7)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                            tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity() + tmp_dist[6].getIntensity();
        double tmp_ave =
            tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
            tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity() + tmp_dist_ave[6].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
        double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
        double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
        double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;
        double tmp_d_7 = tmp_dist[6].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
        double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
        double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
        double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;
        double tmp_a_7 = tmp_dist_ave[6].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist[3].setIntensity(tmp_d_4);
        tmp_dist[4].setIntensity(tmp_d_5);
        tmp_dist[5].setIntensity(tmp_d_6);
        tmp_dist[6].setIntensity(tmp_d_7);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);
        tmp_dist_ave[3].setIntensity(tmp_a_4);
        tmp_dist_ave[4].setIntensity(tmp_a_5);
        tmp_dist_ave[5].setIntensity(tmp_a_6);
        tmp_dist_ave[6].setIntensity(tmp_a_7);

        kld_n7.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(7, tmp_dist, tmp_dist_ave));
      }
    }
    kld_n7.erase(remove(kld_n7.begin(),kld_n7.end(), -9999),kld_n7.end());

    ///N=8

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=8 & tmp_dist_ave.size()>=8)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                            tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity() + tmp_dist[6].getIntensity() + tmp_dist[7].getIntensity();
        double tmp_ave =
            tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
            tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity() + tmp_dist_ave[6].getIntensity() + tmp_dist_ave[7].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
        double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
        double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
        double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;
        double tmp_d_7 = tmp_dist[6].getIntensity() / tmp_digest;
        double tmp_d_8 = tmp_dist[7].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
        double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
        double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
        double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;
        double tmp_a_7 = tmp_dist_ave[6].getIntensity() / tmp_ave;
        double tmp_a_8 = tmp_dist_ave[7].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist[3].setIntensity(tmp_d_4);
        tmp_dist[4].setIntensity(tmp_d_5);
        tmp_dist[5].setIntensity(tmp_d_6);
        tmp_dist[6].setIntensity(tmp_d_7);
        tmp_dist[7].setIntensity(tmp_d_8);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);
        tmp_dist_ave[3].setIntensity(tmp_a_4);
        tmp_dist_ave[4].setIntensity(tmp_a_5);
        tmp_dist_ave[5].setIntensity(tmp_a_6);
        tmp_dist_ave[6].setIntensity(tmp_a_7);
        tmp_dist_ave[7].setIntensity(tmp_a_8);

        kld_n8.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(8, tmp_dist, tmp_dist_ave));
      }
    }
    kld_n8.erase(remove(kld_n8.begin(),kld_n8.end(), -9999),kld_n8.end());

    ///N=9

    for(int i = 0; i < isotope_distributions.size(); ++i)
    {
      IsotopeDistribution tmp_dist = isotope_distributions[i];
      IsotopeDistribution tmp_dist_ave = averagine_distributions[i];

      if(tmp_dist.size()>=9 & tmp_dist_ave.size()>=9)
      {
        double tmp_digest = tmp_dist[0].getIntensity() + tmp_dist[1].getIntensity() + tmp_dist[2].getIntensity() +
                            tmp_dist[3].getIntensity() + tmp_dist[4].getIntensity() + tmp_dist[5].getIntensity() + tmp_dist[6].getIntensity() + tmp_dist[7].getIntensity() + tmp_dist[8].getIntensity();
        double tmp_ave =
            tmp_dist_ave[0].getIntensity() + tmp_dist_ave[1].getIntensity() + tmp_dist_ave[2].getIntensity() +
            tmp_dist_ave[3].getIntensity() + tmp_dist_ave[4].getIntensity() + tmp_dist_ave[5].getIntensity() + tmp_dist_ave[6].getIntensity() + tmp_dist_ave[7].getIntensity()  + tmp_dist_ave[8].getIntensity();


        double tmp_d_1 = tmp_dist[0].getIntensity() / tmp_digest;
        double tmp_d_2 = tmp_dist[1].getIntensity() / tmp_digest;
        double tmp_d_3 = tmp_dist[2].getIntensity() / tmp_digest;
        double tmp_d_4 = tmp_dist[3].getIntensity() / tmp_digest;
        double tmp_d_5 = tmp_dist[4].getIntensity() / tmp_digest;
        double tmp_d_6 = tmp_dist[5].getIntensity() / tmp_digest;
        double tmp_d_7 = tmp_dist[6].getIntensity() / tmp_digest;
        double tmp_d_8 = tmp_dist[7].getIntensity() / tmp_digest;
        double tmp_d_9 = tmp_dist[8].getIntensity() / tmp_digest;

        double tmp_a_1 = tmp_dist_ave[0].getIntensity() / tmp_ave;
        double tmp_a_2 = tmp_dist_ave[1].getIntensity() / tmp_ave;
        double tmp_a_3 = tmp_dist_ave[2].getIntensity() / tmp_ave;
        double tmp_a_4 = tmp_dist_ave[3].getIntensity() / tmp_ave;
        double tmp_a_5 = tmp_dist_ave[4].getIntensity() / tmp_ave;
        double tmp_a_6 = tmp_dist_ave[5].getIntensity() / tmp_ave;
        double tmp_a_7 = tmp_dist_ave[6].getIntensity() / tmp_ave;
        double tmp_a_8 = tmp_dist_ave[7].getIntensity() / tmp_ave;
        double tmp_a_9 = tmp_dist_ave[8].getIntensity() / tmp_ave;

        tmp_dist[0].setIntensity(tmp_d_1);
        tmp_dist[1].setIntensity(tmp_d_2);
        tmp_dist[2].setIntensity(tmp_d_3);
        tmp_dist[3].setIntensity(tmp_d_4);
        tmp_dist[4].setIntensity(tmp_d_5);
        tmp_dist[5].setIntensity(tmp_d_6);
        tmp_dist[6].setIntensity(tmp_d_7);
        tmp_dist[7].setIntensity(tmp_d_8);
        tmp_dist[8].setIntensity(tmp_d_9);
        tmp_dist_ave[0].setIntensity(tmp_a_1);
        tmp_dist_ave[1].setIntensity(tmp_a_2);
        tmp_dist_ave[2].setIntensity(tmp_a_3);
        tmp_dist_ave[3].setIntensity(tmp_a_4);
        tmp_dist_ave[4].setIntensity(tmp_a_5);
        tmp_dist_ave[5].setIntensity(tmp_a_6);
        tmp_dist_ave[6].setIntensity(tmp_a_7);
        tmp_dist_ave[7].setIntensity(tmp_a_8);
        tmp_dist_ave[8].setIntensity(tmp_a_9);


        kld_n9.push_back(DeisotoperRieckert::kullbackLeiblerDivergence(9, tmp_dist, tmp_dist_ave));
      }
    }
    kld_n9.erase(remove(kld_n9.begin(),kld_n9.end(), -9999),kld_n9.end());



    /*
        ///Datenprints für die Ploterstellung, diente der Überprüfung ob es Isotopendistributionen mit den gewünschten KLD Werten gab

        cout << "N2:" <<kld_n2.size() <<endl;

        auto itA2 = max_element(std::begin(kld_n2), std::end(kld_n2));
        cout << *itA2 <<endl;

        auto itB2 = min_element(std::begin(kld_n2), std::end(kld_n2));
        cout << *itB2 <<endl;

        cout << "N3:" <<kld_n3.size() <<endl;

        auto itA3 = max_element(std::begin(kld_n3), std::end(kld_n3));
        cout << *itA3 <<endl;

        auto itB3 = min_element(std::begin(kld_n3), std::end(kld_n3));
        cout << *itB3 <<endl;

        cout << "N4:" <<kld_n4.size() <<endl;

        auto itA4 = max_element(std::begin(kld_n4), std::end(kld_n4));
        cout << *itA4 <<endl;

        auto itB4 = min_element(std::begin(kld_n4), std::end(kld_n4));
        cout << *itB4 <<endl;

        cout << "N5:" <<kld_n5.size() <<endl;

        auto itA5 = max_element(std::begin(kld_n5), std::end(kld_n5));
        cout << *itA5 <<endl;

        auto itB5 = min_element(std::begin(kld_n5), std::end(kld_n5));
        cout << *itB5 <<endl;

        cout << "N6:" <<kld_n6.size() <<endl;

        auto itA6 = max_element(std::begin(kld_n6), std::end(kld_n6));
        cout << *itA6 <<endl;

        auto itB6 = min_element(std::begin(kld_n6), std::end(kld_n6));
        cout << *itB6 <<endl;

        cout << "N7:" <<kld_n7.size() <<endl;

        auto itA7 = max_element(std::begin(kld_n7), std::end(kld_n7));
        cout << *itA7 <<endl;

        auto itB7 = min_element(std::begin(kld_n7), std::end(kld_n7));
        cout << *itB7 <<endl;

        cout << "N8:" <<kld_n8.size() <<endl;

        auto itA8 = max_element(std::begin(kld_n8), std::end(kld_n8));
        cout << *itA8 <<endl;

        auto itB8 = min_element(std::begin(kld_n8), std::end(kld_n8));
        cout << *itB8 <<endl;

        cout << "N9:" <<kld_n9.size() <<endl;

        auto itA9 = max_element(std::begin(kld_n9), std::end(kld_n9));
        cout << *itA9 <<endl;

        auto itB9 = min_element(std::begin(kld_n9), std::end(kld_n9));
        cout << *itB9 <<endl;*/

    ///Schreiben der KL-D Werte in Ausgabedatein zur Ploterstellung
    /*
            ofstream n2stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n2_normalized_realData.txt");
            copy(kld_n2.begin(), kld_n2.end(), std::ostream_iterator<double>(n2stream, "\n"));

            ofstream n3stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n3_normalized_realData.txt");
            copy(kld_n3.begin(), kld_n3.end(), std::ostream_iterator<double>(n3stream, "\n"));

            ofstream n4stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n4_normalized_realData.txt");
            copy(kld_n4.begin(), kld_n4.end(), std::ostream_iterator<double>(n4stream, "\n"));

            ofstream n5stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n5_normalized_realData.txt");
            copy(kld_n5.begin(), kld_n5.end(), std::ostream_iterator<double>(n5stream, "\n"));

            ofstream n6stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n6_normalized_realData.txt");
            copy(kld_n6.begin(), kld_n6.end(), std::ostream_iterator<double>(n6stream, "\n"));

                ofstream n7stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n7_normalized_realData.txt");
                copy(kld_n7.begin(), kld_n7.end(), std::ostream_iterator<double>(n7stream, "\n"));

                ofstream n8stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n8_normalized_realData.txt");
                copy(kld_n8.begin(), kld_n8.end(), std::ostream_iterator<double>(n8stream, "\n"));

                ofstream n9stream("/buffer/ag_bsc/pmsb_2021/Rieckert/data/kldVerify/n9_normalized_realData.txt");
                copy(kld_n9.begin(), kld_n9.end(), std::ostream_iterator<double>(n9stream, "\n"));*/

    cout << "Written Output" << endl;

  }

  ///Berechnung der OpenDeisotoper Zeit

  void DeisotoperRieckert::zeitmessung()
  {
    PeakMap exp_data;
    PeakMap exp_data_picked;
    MzMLFile mzFile;

    mzFile.load("/buffer/ag_bsc/pmsb_2021/Rieckert/data/msdata/MC_bud3_Trypsin_1.mzML",exp_data);

    int peak_count = 0;
    for (auto spec:exp_data)
    {
      for(auto peak:spec)
      {
        peak_count+=1;
      }
    }

    PeakPickerHiRes pp;
    pp.pickExperiment(exp_data,exp_data_picked);

    exp_data_picked.updateRanges();
    vector<MSSpectrum> spectra = exp_data_picked.getSpectra();

    double fragment_tolerance = 5.0;
    bool fragment_unit_ppm = true;
    int min_charge = 1;
    int max_charge = 10;
    bool keep_only_deisotoped = true;
    unsigned int min_isopeaks = 2;
    unsigned int max_isopeaks = 5;
    bool make_single_charged = false;
    bool annotate_charge = true;
    bool annotate_iso_peak_count = true;
    bool use_decreasing_model = true;
    unsigned int start_intensity_check = 1;
    bool add_up_intensity = false;

    auto begin = std::chrono::high_resolution_clock::now();


    for (int i = 0; i < spectra.size(); i++)
    {

      Deisotoper::deisotopeAndSingleCharge(spectra[i],
                                           fragment_tolerance,
                                           fragment_unit_ppm,
                                           min_charge,
                                           max_charge,
                                           keep_only_deisotoped,
                                           min_isopeaks,
                                           max_isopeaks,
                                           make_single_charged,
                                           annotate_charge,
                                           annotate_iso_peak_count,
                                           use_decreasing_model,
                                           start_intensity_check,
                                           add_up_intensity);

    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Für einen Datensatz mit "<<peak_count << " Peaks hat der OpenMS Deisotoper " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "Sekunden gebraucht." << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << "ns" << std::endl;
  }

}