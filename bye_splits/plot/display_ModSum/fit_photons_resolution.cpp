//g++ -o fit_photons_resolution fit_photons_resolution.cpp `root-config --cflags --glibs` -lRooFit -lRooFitCore

//./fit_photons_resolution baseline_photons_-1_5_eta_phi_resolution_12w_9sub.txt baseline -1 PHOTONS CEE_CEH double_sided_crystal_ball
//./fit_photons_resolution 8towers_photons_-1_5_eta_phi_resolution_12w_9sub.txt 8towers_NEW -1 PHOTONS CEE_CEH double_sided_crystal_ball
//./fit_photons_resolution 16towers_photons_-1_5_eta_phi_resolution_12w_9sub.txt 16towers -1 PHOTONS CEE_CEH double_sided_crystal_ball
//./fit_photons_resolution area_overlap_photons_-1_5_eta_phi_resolution_12w_9sub.txt area_overlap -1 PHOTONS CEE_CEH double_sided_crystal_ball

//./fit_photons_resolution 8towers_photons_-1_5_eta_phi_resolution_hadd_123_subw_5x5_OK.txt 8towers -1 PHOTONS CEE_CEH double_sided_crystal_ball


//./fit_photons_resolution filtered_eta_phi_diffs_8towers_CEE_CEH_4k.txt 8towers -1 pions CEE_CEH double_sided_crystal_ball
//baseline_photons_188526_1_eta_phi_resolution_hadd_123_subw_5x5_OK.txt
//baseline_photons_188526_5_eta_phi_resolution_hadd_123_subw_5x5_OK.txt
//16towers_photons_188526_5_eta_phi_resolution_hadd_123_subw_5x5_OK.txt

//baseline_photons_187775_1_eta_phi_resolution_hadd_123_subw_5x5_OK.txt
//double_sided_crystal_ball
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <TCanvas.h>
#include <TH1F.h>
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TLegend.h>

using namespace RooFit;

// Define a custom Double-Sided Crystal Ball function
class RooDoubleSidedCB : public RooAbsPdf {
public:
    RooDoubleSidedCB(const char* name, const char* title,
                     RooAbsReal& x, RooAbsReal& mean, RooAbsReal& sigma,
                     RooAbsReal& alphaL, RooAbsReal& nL,
                     RooAbsReal& alphaR, RooAbsReal& nR)
        : RooAbsPdf(name, title),
          x_("x", "Observable", this, x),
          mean_("mean", "Mean", this, mean),
          sigma_("sigma", "Sigma", this, sigma),
          alphaL_("alphaL", "Alpha Left", this, alphaL),
          nL_("nL", "nL", this, nL),
          alphaR_("alphaR", "Alpha Right", this, alphaR),
          nR_("nR", "nR", this, nR) {}

    RooDoubleSidedCB(const RooDoubleSidedCB& other, const char* name = nullptr)
        : RooAbsPdf(other, name),
          x_("x", this, other.x_),
          mean_("mean", this, other.mean_),
          sigma_("sigma", this, other.sigma_),
          alphaL_("alphaL", this, other.alphaL_),
          nL_("nL", this, other.nL_),
          alphaR_("alphaR", this, other.alphaR_),
          nR_("nR", this, other.nR_) {}

    TObject* clone(const char* newname) const override {
        return new RooDoubleSidedCB(*this, newname);
    }

    inline virtual ~RooDoubleSidedCB() {}

protected:
    RooRealProxy x_;
    RooRealProxy mean_;
    RooRealProxy sigma_;
    RooRealProxy alphaL_;
    RooRealProxy nL_;
    RooRealProxy alphaR_;
    RooRealProxy nR_;

    Double_t evaluate() const override {
        Double_t t = (x_ - mean_) / sigma_;
        if (t < -alphaL_) {
            Double_t a = TMath::Power(nL_ / alphaL_, nL_) * exp(-0.5 * alphaL_ * alphaL_);
            Double_t b = nL_ / alphaL_ - alphaL_;
            return a * TMath::Power(b - t, -nL_);
        } else if (t > alphaR_) {
            Double_t a = TMath::Power(nR_ / alphaR_, nR_) * exp(-0.5 * alphaR_ * alphaR_);
            Double_t b = nR_ / alphaR_ - alphaR_;
            return a * TMath::Power(b + t, -nR_);
        } else {
            return exp(-0.5 * t * t);
        }
    }
};


/*// Updated class for a Double-Sided Crystal Ball
class RooDoubleSidedCB : public RooAbsPdf {
public:
    RooDoubleSidedCB(const char* name, const char* title,
                     const RooAbsReal& x, const RooAbsReal& mean, const RooAbsReal& sigma,
                     const RooAbsReal& alphaL, const RooAbsReal& nL,
                     const RooAbsReal& alphaR, const RooAbsReal& nR)
        : RooAbsPdf(name, title),
          x_("x", "x", this, const_cast<RooAbsReal&>(x)),
          mean_("mean", "mean", this, const_cast<RooAbsReal&>(mean)),
          sigma_("sigma", "sigma", this, const_cast<RooAbsReal&>(sigma)),
          alphaL_("alphaL", "alphaL", this, const_cast<RooAbsReal&>(alphaL)),
          nL_("nL", "nL", this, const_cast<RooAbsReal&>(nL)),
          alphaR_("alphaR", "alphaR", this, const_cast<RooAbsReal&>(alphaR)),
          nR_("nR", "nR", this, const_cast<RooAbsReal&>(nR)) {}

    // Evaluate the PDF
    Double_t evaluate() const override {
        Double_t t = (x_ - mean_) / sigma_;
        if (t < -alphaL_) {
            Double_t a = std::pow(nL_ / alphaL_, nL_) * exp(-0.5 * alphaL_ * alphaL_);
            Double_t b = nL_ / alphaL_ - alphaL_;
            return a * std::pow(b - t, -nL_);
        } else if (t > alphaR_) {
            Double_t a = std::pow(nR_ / alphaR_, nR_) * exp(-0.5 * alphaR_ * alphaR_);
            Double_t b = nR_ / alphaR_ - alphaR_;
            return a * std::pow(b + t, -nR_);
        } else {
            return exp(-0.5 * t * t);
        }
    }

    // Implement the clone method
    TObject* clone(const char* newname) const override {
        return new RooDoubleSidedCB(newname, this->GetTitle(), x_.arg(), mean_.arg(), sigma_.arg(),
                                    alphaL_.arg(), nL_.arg(), alphaR_.arg(), nR_.arg());
    }

private:
    RooRealProxy x_;
    RooRealProxy mean_;
    RooRealProxy sigma_;
    RooRealProxy alphaL_;
    RooRealProxy nL_;
    RooRealProxy alphaR_;
    RooRealProxy nR_;
};*/

void plot_eta_phi_resolution_from_file(const std::string& file_path, const std::string& algo, const std::string& event, const std::string& particle, const std::string& subdet, const std::string& fit_type) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return;
    }

    std::vector<double> eta_diffs;
    std::vector<double> phi_diffs;
    std::string line;

    // Skip the header line
    std::getline(file, line);

    // Read the data
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        double eta_diff, phi_diff;
        char comma;
        ss >> eta_diff >> comma >> phi_diff;
        eta_diffs.push_back(eta_diff);
        phi_diffs.push_back(phi_diff);
    }

    file.close();

    // Create histograms
    TH1F* h_eta = new TH1F("h_eta", "Eta Differences", 25, *std::min_element(eta_diffs.begin(), eta_diffs.end()), *std::max_element(eta_diffs.begin(), eta_diffs.end()));
    TH1F* h_phi = new TH1F("h_phi", "Phi Differences", 25, *std::min_element(phi_diffs.begin(), phi_diffs.end()), *std::max_element(phi_diffs.begin(), phi_diffs.end()));

    for (double val : eta_diffs) h_eta->Fill(val);
    for (double val : phi_diffs) h_phi->Fill(val);

    // Create canvas
    TCanvas* c1 = new TCanvas("c1", "#eta and #phi resolutions", 1200, 600);
    c1->Divide(2, 1);

    // Fit and draw Eta histogram
    c1->cd(1);
    h_eta->Draw();

    //RooRealVar x("x", "Reco #eta - Gen #eta", *std::min_element(eta_diffs.begin(), eta_diffs.end()), *std::max_element(eta_diffs.begin(), eta_diffs.end()));
    RooRealVar x("x", "Reco #eta - Gen #eta", *std::min_element(eta_diffs.begin(), eta_diffs.end()), *std::max_element(eta_diffs.begin(), eta_diffs.end()));
    RooDataHist eta_data("eta_data", "Eta Data", x, Import(*h_eta));

    RooPlot* frame_eta = x.frame();
    eta_data.plotOn(frame_eta);

    if (fit_type == "gaussian") {
        RooRealVar mean_eta("mean_eta", "mean_eta", 0.0, -0.2, 0.2);
        RooRealVar sigma_eta("sigma_eta", "sigma_eta", 0.065, 0.01, 0.2);
        RooGaussian gauss_eta("gauss_eta", "gaussian PDF", x, mean_eta, sigma_eta);
        gauss_eta.fitTo(eta_data);
        gauss_eta.plotOn(frame_eta);

        double rms = h_eta->GetRMS();
        frame_eta->Draw(); //NB: draw the frame before the legend otherwise the legend is not drawn
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.6, 0.8, Form("#mu = %.5f", mean_eta.getVal()));
        latex.DrawLatex(0.6, 0.75, Form("#sigma = %.5f", sigma_eta.getVal()));
        latex.DrawLatex(0.6, 0.70, Form("RMS = %.5f", rms));

    

    } else if (fit_type == "double_sided_crystal_ball") {
        // Define Crystal Ball parameters
        RooRealVar mean_eta("mean_eta", "mean of CB", 0.0, -0.05, 0.05); //-0.001, 0.01
        RooRealVar sigma_eta("sigma_eta", "sigma of CB", 0.05, 0.001, 0.2);
        RooRealVar alphaL_eta("alphaL_eta", "alphaL of CB", 1.5, 0.1, 10);
        RooRealVar alphaR_eta("alphaR_eta", "alphaR of CB", 1.5, 0.1, 10);
        RooRealVar nL_eta("nL_eta", "nL of CB", 2, 0.1, 10);
        RooRealVar nR_eta("nR_eta", "nR of CB", 2, 0.1, 10);
        // Create Crystal Ball PDF
        RooDoubleSidedCB cb_eta("cb_eta", "Crystal Ball PDF", x, mean_eta, sigma_eta, alphaL_eta, nL_eta, alphaR_eta, nR_eta);

        //RooCBShape cb_left("cb_left", "left Crystal Ball", x, mean, sigma, alphaL, nL);
        //RooCBShape cb_right("cb_right", "right Crystal Ball", x, mean, sigma, alphaR, nR);

       // RooAddPdf cb("cb", "double-sided Crystal Ball", RooArgList(cb_left, cb_right));
        double rms_eta = h_eta->GetRMS();
        cb_eta.fitTo(eta_data);
        cb_eta.plotOn(frame_eta);
        frame_eta->Draw();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.65, 0.8, Form("#mu = %.5f", mean_eta.getVal()));
        latex.DrawLatex(0.65, 0.75, Form("#sigma = %.5f", sigma_eta.getVal()));
        latex.DrawLatex(0.65, 0.70, Form("#alpha_{R} = %.5f", alphaR_eta.getVal()));
        latex.DrawLatex(0.65, 0.65, Form("n_{R} = %.5f", nR_eta.getVal()));
        latex.DrawLatex(0.65, 0.60, Form("#alpha_{L} = %.5f", alphaL_eta.getVal()));
        latex.DrawLatex(0.65, 0.55, Form("n_{L} = %.5f", nL_eta.getVal()));
        latex.DrawLatex(0.65, 0.50, Form("RMS = %.5f", rms_eta));
        
        //RooArgSet *nset = new RooArgSet(x);  // Create a normalization set with 'x'
        //cb.fixCoefNormalization(*nset);  // Fix the coefficient normalization using 'x'

        //RooArgSet *nset = new RooArgSet(x);  // Create a normalization set with 'x'
        //cb.fixCoefNormalization(*nset);  // Fix the coefficient normalization using 'x'
    }

    //frame_eta->Draw();
    

    // Fit and draw Phi histogram
    c1->cd(2);
    h_phi->Draw();

    RooRealVar y("y", "Reco #phi - Gen #phi", *std::min_element(phi_diffs.begin(), phi_diffs.end()), *std::max_element(phi_diffs.begin(), phi_diffs.end()));
    RooDataHist phi_data("phi_data", "Phi Data", y, Import(*h_phi));

    RooPlot* frame_phi = y.frame();
    phi_data.plotOn(frame_phi);

    if (fit_type == "gaussian") {
        RooRealVar mean_phi("mean_phi", "mean_phi", 0.0, -0.2, 0.2);
        RooRealVar sigma_phi("sigma_phi", "sigma_phi", 0.08, 0.01, 0.2);
        RooGaussian gauss_phi("gauss_phi", "gaussian PDF", y, mean_phi, sigma_phi);
        gauss_phi.fitTo(phi_data);
        gauss_phi.plotOn(frame_phi);

        double rms_phi = h_phi->GetRMS();
        frame_phi->Draw();
        TLatex latex_phi;
        latex_phi.SetTextSize(0.04);
        latex_phi.SetNDC();
        latex_phi.DrawLatex(0.6, 0.8, Form("#mu = %.5f", mean_phi.getVal()));
        latex_phi.DrawLatex(0.6, 0.75, Form("#sigma = %.5f", sigma_phi.getVal()));
        latex_phi.DrawLatex(0.6, 0.70, Form("RMS = %.5f", rms_phi));


    } else if (fit_type == "double_sided_crystal_ball") {
        RooRealVar mean_phi("mean_phi", "mean_phi", 0.0, -0.001, 0.05);
        RooRealVar sigma_phi("sigma_phi", "sigma_phi", 0.05, 0.001, 0.5);
        RooRealVar alphaL_phi("alphaL_phi", "alphaL_phi", 1.5, 0, 2);
        RooRealVar alphaR_phi("alphaR_phi", "alphaR_phi", 1.5, 0, 2);
        RooRealVar nL_phi("nL_phi", "nL_phi", 1, 0.1, 10);
        RooRealVar nR_phi("nR_phi", "nR_phi", 1, 0.1, 10);
        RooDoubleSidedCB cb_phi("cb_phi", "Crystal Ball PDF phi", y, mean_phi, sigma_phi, alphaL_phi, nL_phi, alphaR_phi, nR_phi);
        //RooCBShape cb_left_phi("cb_left_phi", "left Crystal Ball", y, mean_phi, sigma_phi, alphaL_phi, nL_phi);
        //RooCBShape cb_right_phi("cb_right_phi", "right Crystal Ball", y, mean_phi, sigma_phi, alphaR_phi, nR_phi);

        //RooAddPdf cb_phi("cb_phi", "double-sided Crystal Ball", RooArgList(cb_left_phi, cb_right_phi));
        
        //RooArgSet *nset = new RooArgSet(x);  // Create a normalization set with 'x'
        //cb_left_phi.fixCoefNormalization(*nset);  // Fix the coefficient normalization using 'x'
        double rms_phi = h_phi->GetRMS();

        cb_phi.fitTo(phi_data);
        cb_phi.plotOn(frame_phi);

        frame_phi->Draw();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.65, 0.8, Form("#mu = %.5f", mean_phi.getVal()));
        latex.DrawLatex(0.65, 0.75, Form("#sigma = %.5f", sigma_phi.getVal()));
        latex.DrawLatex(0.65, 0.70, Form("#alpha_{R} = %.5f", alphaR_phi.getVal()));
        latex.DrawLatex(0.65, 0.65, Form("n_{R} = %.5f", nR_phi.getVal()));
        latex.DrawLatex(0.65, 0.60, Form("#alpha_{L} = %.5f", alphaL_phi.getVal()));
        latex.DrawLatex(0.65, 0.55, Form("n_{L} = %.5f", nL_phi.getVal()));
        latex.DrawLatex(0.65, 0.50, Form("RMS = %.5f", rms_phi));

    }

    

    c1->SaveAs((algo + "_" + particle + "_" + event + "_" + subdet + "_"+ fit_type + "_eta_phi_resolution_ROOT.png").c_str());
}

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <file_path> <algo> <event> <particle> <subdet> <fit_type>" << std::endl;
        return 1;
    }

    std::string file_path = argv[1];
    std::string algo = argv[2];
    std::string event = argv[3];
    std::string particle = argv[4];
    std::string subdet = argv[5];
    std::string fit_type = argv[6];

    if (fit_type != "gaussian" && fit_type != "double_sided_crystal_ball") {
        std::cerr << "Error: fit_type must be either 'gaussian' or 'double_sided_crystal_ball'" << std::endl;
        return 1;
    }

    plot_eta_phi_resolution_from_file(file_path, algo, event, particle, subdet, fit_type);

    return 0;
}
