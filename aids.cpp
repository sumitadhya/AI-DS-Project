#include "Pythia8/Pythia.h"
#include <fstream>
#include <random>
#include <cmath>

using namespace Pythia8;

// Function to add Gaussian noise to simulate detector effects
double addNoise(double value, double relativeError) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(value, relativeError * std::fabs(value));
  return d(gen);
}

// Function to process events and write to CSV
void processEvents(Pythia &pythia, std::ofstream &file, int nEvents, int eventClass) {
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    // Event-level aggregated features
    int totalMultiplicity = 0, chargedMultiplicity = 0;
    double sum_pT = 0, sum_pT2 = 0;
    double sum_eta = 0, sum_eta2 = 0;
    double sum_phi = 0, sum_phi2 = 0;
    double totalEnergy = 0;

    for (int i = 0; i < pythia.event.size(); ++i) {
      const auto &particle = pythia.event[i];
      if (!particle.isFinal()) continue;

      // Add noise to simulate detector effects
      double pT = addNoise(particle.pT(), 0.05);  // 5% noise
      double eta = addNoise(particle.eta(), 0.03);
      double phi = addNoise(particle.phi(), 0.02);
      double energy = addNoise(particle.e(), 0.05);

      totalMultiplicity++;
      if (particle.isCharged()) chargedMultiplicity++;

      sum_pT += pT;
      sum_pT2 += pT * pT;
      sum_eta += eta;
      sum_eta2 += eta * eta;
      sum_phi += phi;
      sum_phi2 += phi * phi;
      totalEnergy += energy;
    }

    // Calculate mean and variance
    double mean_pT = sum_pT / totalMultiplicity;
    double var_pT = (sum_pT2 / totalMultiplicity) - (mean_pT * mean_pT);
    double mean_eta = sum_eta / totalMultiplicity;
    double var_eta = (sum_eta2 / totalMultiplicity) - (mean_eta * mean_eta);
    double mean_phi = sum_phi / totalMultiplicity;
    double var_phi = (sum_phi2 / totalMultiplicity) - (mean_phi * mean_phi);

    // Write aggregated features to the CSV
    file << eventClass << "," << totalMultiplicity << "," << chargedMultiplicity << ","
         << mean_pT << "," << var_pT << "," << mean_eta << "," << var_eta << ","
         << mean_phi << "," << var_phi << "," << totalEnergy << "\n";
  }
}

int main() {
  // Initialize Pythia
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000."); // LHC energy at 13 TeV

  // Open CSV files for output
  std::ofstream nonDiffFile("non_diffractive_events.csv");
  std::ofstream diffFile("diffractive_events.csv");

  // Write headers
  nonDiffFile << "Class,TotalMultiplicity,ChargedMultiplicity,Mean_pT,Var_pT,Mean_Eta,Var_Eta,Mean_Phi,Var_Phi,TotalEnergy\n";
  diffFile << "Class,TotalMultiplicity,ChargedMultiplicity,Mean_pT,Var_pT,Mean_Eta,Var_Eta,Mean_Phi,Var_Phi,TotalEnergy\n";

  int nEvents = 1000; // Number of events per class

  // ** Non-Diffractive Events **
  pythia.readString("HardQCD:all = on"); // Enable non-diffractive processes
  pythia.init();
  processEvents(pythia, nonDiffFile, nEvents, 0); // Label: 0 for non-diffractive

  // ** Diffractive Events **
  pythia.readString("HardQCD:all = off"); // Disable non-diffractive processes
  pythia.readString("SoftQCD:singleDiffractive = on");
  pythia.readString("SoftQCD:doubleDiffractive = on");
  pythia.readString("SoftQCD:centralDiffractive = on");
  pythia.init();
  processEvents(pythia, diffFile, nEvents, 1); // Label: 1 for diffractive

  // Close files and print statistics
  nonDiffFile.close();
  diffFile.close();
  pythia.stat();

  return 0;
}
