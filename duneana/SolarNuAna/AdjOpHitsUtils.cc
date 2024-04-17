#include "AdjOpHitsUtils.h"

namespace solar
{
  AdjOpHitsUtils::AdjOpHitsUtils(fhicl::ParameterSet const &p)
      : fOpFlashAlgoTime(p.get<double>("OpFlashAlgoTime")),
        fOpFlashAlgoRad(p.get<double>("OpFlashAlgoRad")),
        fAdjOpFlashMinPECut(p.get<float>("AdjOpFlashMinPECut"))
  {
  }
  void AdjOpHitsUtils::MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, art::Event const &evt){
    // This is the constructor for the OpFlash object
    // OpFlash (double time, double timewidth, double abstime, unsigned int frame, std::vector< double > PEperOpDet, bool InBeamFrame=0, int OnBeamTime=0, double FastToTotal=1, double yCenter=0, double yWidth=0, double zCenter=0, double zWidth=0, std::vector< double > WireCenters=std::vector< double >(0), std::vector< double > WireWidths=std::vector< double >(0))
    auto assns = std::make_unique<art::Assns<recob::OpFlash, recob::OpHit>>(); 
    
    for (std::vector<art::Ptr<recob::OpHit>> Cluster : Clusters)
    {
      std::vector<art::Ptr<recob::OpHit>> ClusterCopy = Cluster;
      std::sort(ClusterCopy.begin(), ClusterCopy.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b) { return a->PeakTime() < b->PeakTime(); });
      int NHit = 0;
      double Time = ClusterCopy[0]->PeakTime();
      double TimeWidth = 0;
      double TimeSum = 0;
      double PE = 0;
      double MaxPE = 0;
      std::vector<double> PEperOpDet = {};
      double FastToTotal = 1;
      double X = 0;
      double Y = 0;
      double Z = 0;
      double YWidth = 0;
      double ZWidth = 0;
      double XSum = 0;
      double YSum = 0;
      double ZSum = 0;
      
      for (art::Ptr<recob::OpHit> PDSHit : ClusterCopy)
      {
        NHit++;
        TimeSum += PDSHit->PeakTime();
        PE += PDSHit->PE();
        if (PDSHit->PE() > MaxPE)
          MaxPE = PDSHit->PE();
        PEperOpDet.push_back(PDSHit->PE());
        auto OpHitXYZ = geo->OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        XSum += OpHitXYZ.X();
        YSum += OpHitXYZ.Y();
        ZSum += OpHitXYZ.Z();
      }
      Time = TimeSum / ClusterCopy.size();
      X = XSum / ClusterCopy.size();
      Y = YSum / ClusterCopy.size();
      Z = ZSum / ClusterCopy.size();
      for (art::Ptr<recob::OpHit> PDSHit : ClusterCopy)
      {
        auto OpHitXYZ = geo->OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        TimeWidth += (PDSHit->PeakTime() - Time) * (PDSHit->PeakTime() - Time);
        YWidth += (OpHitXYZ.Y() - Y) * (OpHitXYZ.Y() - Y);
        ZWidth += (OpHitXYZ.Z() - Z) * (OpHitXYZ.Z() - Z);
      }
      TimeWidth = sqrt(TimeWidth / ClusterCopy.size());
      YWidth = sqrt(YWidth / ClusterCopy.size());
      ZWidth = sqrt(ZWidth / ClusterCopy.size());
      // Compute FastToTotal according to the #PEs arriving within the first 10% of the time window wrt the total #PEs
      for (art::Ptr<recob::OpHit> PDSHit : ClusterCopy)
      {
        if (PDSHit->PeakTime() < Time + TimeWidth / 10)
          FastToTotal += PDSHit->PE();
      }
      FastToTotal /= PE;
      FlashVec.push_back(FlashInfo{NHit, Time, TimeWidth, PE, MaxPE, PEperOpDet, FastToTotal, X, Y, Z, YWidth, ZWidth});
      // Generate the flash object according to the constructor
      // FlashInfo{Time, TimeWidth, PE, PEperOpDet, FastToTotal, Y, Z, YWidth, ZWidth};
      // recob::OpFlash Flash(Time, TimeWidth, Time, 0, PEperOpDet, 0, 0, FastToTotal, Y, YWidt.h, Z, ZWidth);
      // auto const FlashPtrMaker = art::PtrMaker<recob::OpFlash>(evt);
      // art::Ptr<recob::OpFlash> FlashPtr = FlashPtrMaker(FlashVec.size()-1);
      // for (art::Ptr<recob::OpHit> PDSHit : Cluster)
      // {
      //   assns->addSingle(FlashPtr, PDSHit);
      // }
    }
    // evt.put(std::move(assns));
    return;
  }
  void AdjOpHitsUtils::CalcAdjOpHits(std::vector<art::Ptr<recob::OpHit>> MyVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, bool HeavDebug)
  /*
  Find adjacent hits in time and space:
  - MyVec is the vector of hits to be clustered
  - Clusters is the vector of clusters
  - MyHist is the histogram to be filled with the number of hits in each cluster
  - ADCIntHist is the histogram to be filled with the ADC integral of each cluster
  - HeavDebug is a boolean to turn on/off debugging statements
  */
  {
    const double TimeRange = fOpFlashAlgoTime; // Time in ns
    const double RadRange = fOpFlashAlgoRad; // Range in cm
    unsigned int FilledHits = 0;
    unsigned int NumOriHits = MyVec.size();

    while (NumOriHits != FilledHits)
    {
      if (HeavDebug)
        std::cerr << "\nStart of my while loop" << std::endl;
      std::vector<art::Ptr<recob::OpHit>> AdjHitVec;
      AdjHitVec.push_back(MyVec[0]);
      MyVec.erase(MyVec.begin() + 0);
      int LastSize = 0;
      int NewSize = AdjHitVec.size();

      while (LastSize != NewSize)
      {
        std::vector<int> AddNow;
        for (size_t aL = 0; aL < AdjHitVec.size(); ++aL)
        {
          for (size_t nL = 0; nL < MyVec.size(); ++nL)
          {
            if (HeavDebug)
            {
              std::cerr << "\t\tLooping though AdjVec " << aL << " and  MyVec " << nL
                        << " AdjHitVec - " << AdjHitVec[aL]->OpChannel() << " & " << AdjHitVec[aL]->PeakTime()
                        << " MVec - " << MyVec[nL]->OpChannel() << " & " << MyVec[nL]->PeakTime()
                        // << " OpChannel " << abs((int)AdjHitVec[aL]->OpChannel() - (int)MyVec[nL]->OpChannel()) << " bool " << (bool)(abs((int)AdjHitVec[aL]->OpChannel() - (int)MyVec[nL]->OpChannel()) <= ChanRange)
                        << " Time " << abs(AdjHitVec[aL]->PeakTime() - MyVec[nL]->PeakTime()) << " bool " << (bool)(abs((double)AdjHitVec[aL]->PeakTime() - (double)MyVec[nL]->PeakTime()) <= TimeRange)
                        << std::endl;
            }
            auto MyHitXYZ = geo->OpDetGeoFromOpChannel((int)MyVec[nL]->OpChannel()).GetCenter();
            auto AdjHitXYZ = geo->OpDetGeoFromOpChannel((int)AdjHitVec[aL]->OpChannel()).GetCenter();
            if ((AdjHitXYZ - MyHitXYZ).R() <= RadRange &&
              abs((double)AdjHitVec[aL]->PeakTime() - (double)MyVec[nL]->PeakTime()) <= TimeRange)
            {

              if (HeavDebug)
                std::cerr << "\t\t\tFound a new thing!!!" << std::endl;
              // --- Check that this element isn't already in AddNow.
              bool AlreadyPres = false;

              for (size_t zz = 0; zz < AddNow.size(); ++zz)
              {
                if (AddNow[zz] == (int)nL)
                  AlreadyPres = true;
              }

              if (!AlreadyPres)
                AddNow.push_back(nL);
            } // If this PDSHit is within the window around one of my other hits.
          }   // Loop through my vector of collection plane hits.
        }     // Loop through AdjHitVec

        // --- Now loop through AddNow and remove from Marley whilst adding to AdjHitVec
        std::sort(AddNow.begin(), AddNow.end());
        for (size_t aa = 0; aa < AddNow.size(); ++aa)
        {
          if (HeavDebug)
          {
            std::cerr << "\tRemoving element " << AddNow.size() - 1 - aa << " from MyVec ===> "
                      << MyVec[AddNow[AddNow.size() - 1 - aa]]->OpChannel() << " & " << MyVec[AddNow[AddNow.size() - 1 - aa]]->PeakTime()
                      << std::endl;
          }

          AdjHitVec.push_back(MyVec[AddNow[AddNow.size() - 1 - aa]]);
          MyVec.erase(MyVec.begin() + AddNow[AddNow.size() - 1 - aa]); // This line creates segmentation fault
                                                                      // std::cout << "Erase works" << std::endl;
        }

        LastSize = NewSize;
        NewSize = AdjHitVec.size();
        if (HeavDebug)
        {
          std::cerr << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
                    << "\nLets see what is in AdjHitVec...." << std::endl;
          for (size_t aL = 0; aL < AdjHitVec.size(); ++aL)
          {
            std::cout << "\tElement " << aL << " is ===> " << AdjHitVec[aL]->OpChannel() << " & " << AdjHitVec[aL]->PeakTime() << std::endl;
          }
        }
      } // while ( LastSize != NewSize )

      int NumAdjColHits = AdjHitVec.size();
      float SummedADCInt = 0;
      for (art::Ptr<recob::OpHit> PDSHit : AdjHitVec)
        SummedADCInt += PDSHit->PE();

      if (HeavDebug)
        std::cerr << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits." << std::endl;

      FilledHits += NumAdjColHits;

      if (AdjHitVec.size() > 0)
        Clusters.push_back(AdjHitVec);
    }

    if (HeavDebug)
    {
      std::vector<double> avgChannel;
      std::vector<double> avgTick;
      std::vector<double> summedADCInt;

      for (std::vector<art::Ptr<recob::OpHit>> hits : Clusters)
      {
        double adcInt = 0;
        double channel = 0;
        double tick = 0;

        for (art::Ptr<recob::OpHit> PDSHit : hits)
        {
          tick += PDSHit->PE() * PDSHit->PeakTime();
          channel += PDSHit->PE() * PDSHit->OpChannel();
          adcInt += PDSHit->PE();
        }
        tick /= adcInt;
        channel /= adcInt;
        summedADCInt.push_back(adcInt);
        avgTick.push_back(tick);
        avgChannel.push_back(channel);
      }

      for (int i = 0; i < int(avgTick.size() - 1); i++)
      {
        for (int j = i + 1; j < int(avgTick.size()); j++)
        {
          std::cout << avgChannel[i] << " " << avgChannel[j] << "  " << std::abs(avgChannel[i] - avgChannel[j]) << std::endl;
          std::cout << avgTick[i] << " " << avgTick[j] << "  " << std::abs(avgTick[i] - avgTick[j]) << std::endl;
          std::cout << summedADCInt[i] << " " << summedADCInt[j] << std::endl;
        }
      }
    }
    return;
  }
} // namespace solar