#include "Pythia8/Pythia.h"
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> // for formatting output

using namespace Pythia8;

// Function to calculate rapidity
double calculateRapidity(double energy, double pz) {
  return 0.5 * std::log((energy + pz) / (energy - pz));
}

// Function to process events and write particle data to CSV
void processEvents(Pythia &pythia, std::ofstream &file, int nEvents, int eventClass) {
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    for (int i = 0; i < pythia.event.size(); ++i) {
      const auto &particle = pythia.event[i];
      if (!particle.isFinal()) continue; // Only consider final-state particles

      // Extract particle properties
      double pT = particle.pT();
      double eta = particle.eta();
      double phi = particle.phi();
      double energy = particle.e();
      double pz = particle.pz();
      double rapidity = calculateRapidity(energy, pz);

      // Write particle data to the CSV file
      file << iEvent << ","      // Event number
           << particle.id() << "," // Particle ID (PDG code)
           << particle.isCharged() << "," // Charged (1) or neutral (0)
           << pT << ","          // Transverse momentum
           << eta << ","         // Pseudorapidity
           << phi << ","         // Azimuthal angle
           << rapidity << ","    // Rapidity
           << energy << ","      // Energy
           << eventClass << "\n"; // Event class (0: non-diffractive, 1: diffractive)
    }
  }
}

int main() {
  // Initialize Pythia
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000."); // LHC energy at 13 TeV
  pythia.readString("Beams:idA = 2212");   // Proton beam 1
  pythia.readString("Beams:idB = 2212");   // Proton beam 2

  // Open CSV files for output
  std::ofstream nonDiffFile("non_diffractive_particles.csv");
  std::ofstream diffFile("diffractive_particles.csv");

  // Write headers
  nonDiffFile << "Event,ParticleID,Charged,pT,Eta,Phi,Rapidity,Energy,Class\n";
  diffFile << "Event,ParticleID,Charged,pT,Eta,Phi,Rapidity,Energy,Class\n";

  int nEvents = 5000; // Number of events per class

  // ** Non-Diffractive Events **
  pythia.readString("HardQCD:all = on"); // Enable non-diffractive processes
  pythia.init();
  processEvents(pythia, nonDiffFile, nEvents, 0); // Label: 0 for non-diffractive
  pythia.stat(); // Print statistics for non-diffractive events
  pythia.settings.resetAll(); // Reset settings before switching to diffractive events

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
