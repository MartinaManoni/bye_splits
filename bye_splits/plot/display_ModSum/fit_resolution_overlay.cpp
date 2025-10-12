// g++ -o fit_resolution_overlay fit_resolution_overlay.cpp `root-config --cflags --glibs` -lRooFit -lRooFitCore -lASImage
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <TCanvas.h>
#include <TH1F.h>
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
#include <TPaveText.h>

// CMS style macros
#include "tdrstyle.C"
#include "CMS_lumi.C"

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

void plot_eta_phi_resolution_from_file(
    const std::string& file_path1,
    const std::string& file_path2,
    const std::string& algo,
    const std::string& event,
    const std::string& particle,
    const std::string& subdet,
    const std::string& fit_type) {

    // --- CMS style ---
    setTDRStyle();
    int iPeriod = 0;   // 13 TeV
    int iPos    = 0;   // top right

    // --- Read both files ---
    std::vector<double> eta_diffs, phi_diffs;
    std::vector<double> eta_diffs2, phi_diffs2;
    std::string line;

    auto read_file = [](const std::string& path, std::vector<double>& eta, std::vector<double>& phi){
        std::ifstream file(path);
        if (!file.is_open()) { 
            std::cerr << "Error opening file " << path << std::endl; 
            return; 
        }

        std::string line; // <-- declare here
        std::getline(file, line); // skip header

        double e, p; char comma;
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            ss >> e >> comma >> p;
            eta.push_back(e);
            phi.push_back(p);
        }
        file.close();
    };

    read_file(file_path1, eta_diffs, phi_diffs);
    read_file(file_path2, eta_diffs2, phi_diffs2);


    // --- Combined ranges ---
    double x_min = std::min(*std::min_element(eta_diffs.begin(), eta_diffs.end()),
                            *std::min_element(eta_diffs2.begin(), eta_diffs2.end()));
    double x_max = std::max(*std::max_element(eta_diffs.begin(), eta_diffs.end()),
                            *std::max_element(eta_diffs2.begin(), eta_diffs2.end()));
    double y_min = std::min(*std::min_element(phi_diffs.begin(), phi_diffs.end()),
                            *std::min_element(phi_diffs2.begin(), phi_diffs2.end()));
    double y_max = std::max(*std::max_element(phi_diffs.begin(), phi_diffs.end()),
                            *std::max_element(phi_diffs2.begin(), phi_diffs2.end()));

    // --- Create histograms ---
    TH1F* h_eta  = new TH1F("h_eta", ";Reco #eta - Gen #eta;Events", 25, x_min, x_max);
    TH1F* h_eta2 = new TH1F("h_eta2", ";Reco #eta - Gen #eta;Events", 25, x_min, x_max);
    TH1F* h_phi  = new TH1F("h_phi", ";Reco #phi - Gen #phi;Events", 25, y_min, y_max);
    TH1F* h_phi2 = new TH1F("h_phi2", ";Reco #phi - Gen #phi;Events", 25, y_min, y_max);


    for (double v : eta_diffs) h_eta->Fill(v);
    for (double v : eta_diffs2) h_eta2->Fill(v);
    for (double v : phi_diffs) h_phi->Fill(v);
    for (double v : phi_diffs2) h_phi2->Fill(v);

    // --- Eta ---
    TCanvas* c_eta = new TCanvas("c_eta", "#eta Resolution", 800, 600);
    //c_eta->SetLeftMargin(0.15);   // default ~0.13, increase to add more space for y-axis title
    //c_eta->SetRightMargin(0.05);  // default ~0.05
    c_eta->SetTopMargin(0.1);    // default ~0.08
    //c_eta->SetBottomMargin(0.15); // default ~0.13, increase to add space for x-axis title  
    RooRealVar x("x", "Reco #eta - Gen #eta", x_min, x_max);

    RooDataHist eta_data("eta_data", "Eta Data", x, Import(*h_eta));
    RooDataHist eta_data2("eta_data2", "Eta Data2", x, Import(*h_eta2));

    RooPlot* frame_eta = x.frame();
    //eta_data.plotOn(frame_eta);
    //eta_data2.plotOn(frame_eta, LineColor(kRed));

    eta_data.plotOn(frame_eta, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(20));
    eta_data2.plotOn(frame_eta, MarkerColor(kRed), LineColor(kRed), MarkerStyle(20));

    double binWidth_eta = h_eta->GetBinWidth(1);
    TString ytitle_eta = Form("Events / (%.3f)", binWidth_eta); // keep 3 decimals
    frame_eta->GetYaxis()->SetTitle(ytitle_eta);
    frame_eta->GetYaxis()->SetTitleOffset(1.00);
    frame_eta->GetXaxis()->SetLabelSize(0.05);
    frame_eta->GetYaxis()->SetLabelSize(0.05);

    // TPaveText labels
    TPaveText* label1 = new TPaveText(0.20, 0.75, 0.48, 0.88, "NDC");
    TPaveText* label2 = new TPaveText(0.65, 0.75, 0.90, 0.88, "NDC");

    TPaveText* label3 = new TPaveText(0.20, 0.75, 0.48, 0.88, "NDC");
    TPaveText* label4 = new TPaveText(0.65, 0.75, 0.90, 0.88, "NDC");

    label1->SetFillColor(0); label1->SetFillStyle(0); label1->SetLineColor(0); label1->SetShadowColor(0);
    label2->SetFillColor(0); label2->SetFillStyle(0); label2->SetLineColor(0); label2->SetShadowColor(0);

    label3->SetFillColor(0); label3->SetFillStyle(0); label3->SetLineColor(0); label3->SetShadowColor(0);
    label4->SetFillColor(0); label4->SetFillStyle(0); label4->SetLineColor(0); label4->SetShadowColor(0);

    label1->SetFillColor(0); label2->SetFillColor(0);
    label1->SetTextAlign(12); label2->SetTextAlign(12);
    label1->SetTextSize(0.04); label2->SetTextSize(0.04);
    label1->SetTextFont(42);  label2->SetTextFont(42);

    label3->SetFillColor(0); label4->SetFillColor(0);
    label3->SetTextAlign(12); label4->SetTextAlign(12);
    label3->SetTextSize(0.04); label4->SetTextSize(0.04);
    label3->SetTextFont(42);  label4->SetTextFont(42);
    

    if (fit_type == "gaussian") {
        RooRealVar mean_eta("mean_eta", "mean_eta", 0.0, -0.2, 0.2);
        RooRealVar sigma_eta("sigma_eta", "sigma_eta", 0.065, 0.01, 0.2);
        RooGaussian gauss_eta("gauss_eta", "gaussian PDF", x, mean_eta, sigma_eta);

        gauss_eta.fitTo(eta_data, Range(x_min, x_max), Extended(kFALSE));
        gauss_eta.plotOn(frame_eta);
        gauss_eta.fitTo(eta_data2, Range(x_min, x_max), Extended(kFALSE));
        gauss_eta.plotOn(frame_eta, LineColor(kRed));

    } else {
        RooRealVar mean_eta("mean_eta", "mean of CB", 0.0, -0.05, 0.05);
        RooRealVar sigma_eta("sigma_eta", "sigma of CB", 0.05, 0.001, 0.2);
        RooRealVar alphaL_eta("alphaL_eta", "alphaL of CB", 1.5, 0.1, 10);
        RooRealVar alphaR_eta("alphaR_eta", "alphaR of CB", 1.5, 0.1, 10);
        RooRealVar nL_eta("nL_eta", "nL of CB", 2, 0.1, 10);
        RooRealVar nR_eta("nR_eta", "nR of CB", 2, 0.1, 10);
        RooDoubleSidedCB cb_eta("cb_eta", "Crystal Ball PDF", x,
                                mean_eta, sigma_eta,
                                alphaL_eta, nL_eta, alphaR_eta, nR_eta);

        cb_eta.fitTo(eta_data, Range(x_min, x_max), Extended(kFALSE));

        // Extract parameters for dataset 1
        double mean1 = mean_eta.getVal();
        double sigma1 = sigma_eta.getVal();
        double alphaL1 = alphaL_eta.getVal();
        double alphaR1 = alphaR_eta.getVal();
        double nL1 = nL_eta.getVal();
        double nR1 = nR_eta.getVal();
        //double chi2_1 = frame_eta->chiSquare();  

        cb_eta.plotOn(frame_eta, Name("fit1"));
        double chi2_1 = frame_eta->chiSquare("fit1", "eta_data",6);
        cb_eta.fitTo(eta_data2, Range(x_min, x_max), Extended(kFALSE));

        // Extract parameters for dataset 2
        double mean2 = mean_eta.getVal();
        double sigma2 = sigma_eta.getVal();
        double alphaL2 = alphaL_eta.getVal();
        double alphaR2 = alphaR_eta.getVal();
        double nL2 = nL_eta.getVal();
        double nR2 = nR_eta.getVal();
        //double chi2_2 = frame_eta->chiSquare();  

        cb_eta.plotOn(frame_eta, LineColor(kRed), Name("fit2"));
        double chi2_2 = frame_eta->chiSquare("fit2", "eta_data2",6);

        // Add labels line by line
        label1->SetTextColor(kBlue);
        TText* t1 = label1->AddText("No module splitting");
        t1->SetTextFont(62); 
        //label1->AddText(Form("#mu = %.3f", mean1));
        //label1->AddText(Form("#sigma = %.3f", sigma1));
        //label1->AddText(Form("#alpha_{L} = %.3f", alphaL1));
        //label1->AddText(Form("n_{L} = %.3f", nL1));
        //label1->AddText(Form("#alpha_{R} = %.3f", alphaR1));
        //label1->AddText(Form("n_{R} = %.3f", nR1));

        //label1->AddText("#sigma_{eff} =  0.017"); //Baseline PIONS
        label1->AddText("#sigma_{eff} =  0.025"); //Baseline Jets

        //label1->AddText(Form("#chi^{2}/NDF = %.3f", chi2_1));

        label2->SetTextColor(kRed);
        TText* t2 = label2->AddText("1/16 module splitting");
        t2->SetTextFont(62); 
        //label2->AddText(Form("#mu = %.3f", mean2));
        //label2->AddText(Form("#sigma = %.3f", sigma2));
        //label2->AddText(Form("#alpha_{L} = %.3f", alphaL2));
        //label2->AddText(Form("n_{L} = %.3f", nL2));
        //label2->AddText(Form("#alpha_{R} = %.3f", alphaR2));
        //label2->AddText(Form("n_{R} = %.3f", nR2));


        //label2->AddText("#sigma_{eff} =  0.017"); //1/16 PIONS
        label2->AddText("#sigma_{eff} =  0.025"); //1/16 Jets



        //label2->AddText(Form("#chi^{2}/NDF = %.3f", chi2_2));
        
    }

    frame_eta->Draw();
    label1->Draw();
    label2->Draw();
    // Lower-right text, in NDC coordinates (0–1)
    double xNDC0 = 0.21;
    double yNDC0 = 0.60;
    TLatex* t30 = new TLatex(xNDC0, yNDC0, "Jets PU=0");
    t30->SetNDC();              // use normalized coordinates
    t30->SetTextColor(kBlack);
    t30->SetTextFont(42);
    t30->SetTextAlign(11);      // right-aligned horizontally
    t30->Draw();

    CMS_lumi(c_eta, 0, 0);
    c_eta->SaveAs((algo + "_" + particle + "_" + event + "_" + subdet + "_" + fit_type + "_eta_overlay.png").c_str());
    c_eta->SaveAs((algo + "_" + particle + "_" + event + "_" + subdet + "_" + fit_type + "_eta_overlay.pdf").c_str());

    // --- Phi ---
    TCanvas* c_phi = new TCanvas("c_phi", "#phi Resolution", 800, 600);
    c_phi->SetTopMargin(0.1);
    RooRealVar y("y", "Reco #phi - Gen #phi", y_min, y_max);
    RooDataHist phi_data("phi_data", "Phi Data", y, Import(*h_phi));
    RooDataHist phi_data2("phi_data2", "Phi Data2", y, Import(*h_phi2));

    RooPlot* frame_phi = y.frame();
    phi_data.plotOn(frame_phi, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(20));
    phi_data2.plotOn(frame_phi, MarkerColor(kRed), LineColor(kRed), MarkerStyle(20));

    double binWidth_phi = h_phi->GetBinWidth(1);
    TString ytitle_phi = Form("Events / (%.3f)", binWidth_phi); // keep 3 decimals
    frame_phi->GetYaxis()->SetTitle(ytitle_phi);
    frame_phi->GetYaxis()->SetTitleOffset(1.00);
    frame_phi->GetXaxis()->SetLabelSize(0.05);
    frame_phi->GetYaxis()->SetLabelSize(0.05);

    if (fit_type == "gaussian") {
        RooRealVar mean_phi("mean_phi", "mean_phi", 0.0, -0.2, 0.2);
        RooRealVar sigma_phi("sigma_phi", "sigma_phi", 0.08, 0.01, 0.2);
        RooGaussian gauss_phi("gauss_phi", "gaussian PDF", y, mean_phi, sigma_phi);

        gauss_phi.fitTo(phi_data, Range(y_min, y_max), Extended(kFALSE));
        gauss_phi.plotOn(frame_phi);
        gauss_phi.fitTo(phi_data2, Range(y_min, y_max), Extended(kFALSE));
        gauss_phi.plotOn(frame_phi, LineColor(kRed));

    } else {
        RooRealVar mean_phi("mean_phi", "mean_phi", 0.0, -0.001, 0.05);
        RooRealVar sigma_phi("sigma_phi", "sigma_phi", 0.05, 0.001, 0.5);
        RooRealVar alphaL_phi("alphaL_phi", "alphaL_phi", 1.5, 0, 2);
        RooRealVar alphaR_phi("alphaR_phi", "alphaR_phi", 1.5, 0, 2);
        RooRealVar nL_phi("nL_phi", "nL_phi", 1, 0.1, 10);
        RooRealVar nR_phi("nR_phi", "nR_phi", 1, 0.1, 10);
        RooDoubleSidedCB cb_phi("cb_phi", "Crystal Ball PDF phi", y,
                                mean_phi, sigma_phi,
                                alphaL_phi, nL_phi, alphaR_phi, nR_phi);

        cb_phi.fitTo(phi_data, Range(y_min, y_max), Extended(kFALSE));

        // Extract parameters for dataset 1
        double mean1 = mean_phi.getVal();
        double sigma1 = sigma_phi.getVal();
        double alphaL1 = alphaL_phi.getVal();
        double alphaR1 = alphaR_phi.getVal();
        double nL1 = nL_phi.getVal();
        double nR1 = nR_phi.getVal();
        //double chi2_phi = frame_phi->chiSquare();
        //double chi2_1 = frame_phi->chiSquare(); 

        cb_phi.plotOn(frame_phi);

        cb_phi.fitTo(phi_data2, Range(y_min, y_max), Extended(kFALSE));
        

        // Extract parameters for dataset 2
        double mean2 = mean_phi.getVal();
        double sigma2 = sigma_phi.getVal();
        double alphaL2 = alphaL_phi.getVal();
        double alphaR2 = alphaR_phi.getVal();
        double nL2 = nL_phi.getVal();
        double nR2 = nR_phi.getVal();

        cb_phi.plotOn(frame_phi, LineColor(kRed));


        // Add labels line by line
        label3->SetTextColor(kBlue);
        TText* t1 = label3->AddText("No module splitting");
        t1->SetTextFont(62); 
        //label3->AddText(Form("#mu = %.3f", mean1));
        //label3->AddText(Form("#sigma = %.3f", sigma1));
        //label3->AddText(Form("#alpha_{L} = %.3f", alphaL1));
        //label3->AddText(Form("n_{L} = %.3f", nL1));
        //label3->AddText(Form("#alpha_{R} = %.3f", alphaR1));
        //label3->AddText(Form("n_{R} = %.3f", nR1));

        //label3->AddText("#sigma_{eff} =  0.021"); //Baseline PIONS
        label3->AddText("#sigma_{eff} =  0.030"); //Baseline Jets

        //label3->AddText(Form("#chi^{2}/NDF = %.3f", chi2_1));

        label4->SetTextColor(kRed);
        TText* t2 = label4->AddText("1/16 module splitting");
        t2->SetTextFont(62); 
        //label4->AddText(Form("#mu = %.3f", mean2));
        //label4->AddText(Form("#sigma = %.3f", sigma2));
        //label4->AddText(Form("#alpha_{L} = %.3f", alphaL2));
        //label4->AddText(Form("n_{L} = %.3f", nL2));
        //label4->AddText(Form("#alpha_{R} = %.3f", alphaR2));
        //label4->AddText(Form("n_{R} = %.3f", nR2));

        //label4->AddText("#sigma_{eff} =  0.019"); //1/16 PIONS
        label4->AddText("#sigma_{eff} =  0.027"); //1/16 Jets

        //label4->AddText(Form("#chi^{2}/NDF = %.3f", chi2_2));
    }

    frame_phi->Draw();
    label3->Draw();
    label4->Draw();


    // Lower-right text, in NDC coordinates (0–1)
    double xNDC = 0.21;
    double yNDC = 0.60;
    TLatex* t3 = new TLatex(xNDC, yNDC, "Jets PU=0");
    t3->SetNDC();              // use normalized coordinates
    t3->SetTextColor(kBlack);
    t3->SetTextFont(42);
    t3->SetTextAlign(11);      // right-aligned horizontally
    t3->Draw();

    CMS_lumi(c_phi, 0, 0);
    c_phi->SaveAs((algo + "_" + particle + "_" + event + "_" + subdet + "_" + fit_type + "_phi_overlay.png").c_str());
    c_phi->SaveAs((algo + "_" + particle + "_" + event + "_" + subdet + "_" + fit_type + "_phi_overlay.pdf").c_str());
}



int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] 
                  << " <file_path1> <file_path2> <algo> <event> <particle> <subdet> <fit_type>" 
                  << std::endl;
        return 1;
    }

    std::string file_path1 = argv[1];   // first file
    std::string file_path2 = argv[2];   // second file
    std::string algo       = argv[3];
    std::string event      = argv[4];
    std::string particle   = argv[5];
    std::string subdet     = argv[6];
    std::string fit_type   = argv[7];

    if (fit_type != "gaussian" && fit_type != "double_sided_crystal_ball") {
        std::cerr << "Error: fit_type must be either 'gaussian' or 'double_sided_crystal_ball'" << std::endl;
        return 1;
    }

    plot_eta_phi_resolution_from_file(file_path1, file_path2, algo, event, particle, subdet, fit_type);

    return 0;
}
