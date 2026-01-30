Double_t GetEnergy(Double_t P, Double_t m);     // E
Double_t GetRapidity(Double_t E, Double_t Pz);  // y =Atanh (beta)
Double_t GetTransMom(Double_t Px, Double_t Py); // pT=sqrt(px^2+py^2)
Double_t GetTotMom(Double_t Px, Double_t Py,
                   Double_t Pz);                // p=sqrt(px^2+py^2+pz^2)
Double_t GetTransMass(Double_t m, Double_t pT); // mT =  sqrt(pT^2+m^2)
Double_t GetPhi(Double_t fpx, Double_t fpy);    // phi = ATan(py/px)
Double_t GetEta(Double_t theta);
Double_t GetTheta(Double_t Pz, Double_t P);
Double_t GetV2(Double_t phi);
Double_t GetV3(Double_t phi);

void sethisto(TH2D *Graph_Graph03001)
{
  Graph_Graph03001->SetDirectory(0);
  Graph_Graph03001->SetStats(0);

  Graph_Graph03001->GetXaxis()->SetLabelFont(42);
  Graph_Graph03001->GetXaxis()->SetLabelOffset(0.004);
  Graph_Graph03001->GetXaxis()->SetLabelSize(0.045);
  Graph_Graph03001->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph03001->GetXaxis()->SetTickLength(0.03);
  Graph_Graph03001->GetXaxis()->SetTitleFont(42);
  Graph_Graph03001->GetXaxis()->SetNdivisions(510, kTRUE);
  Graph_Graph03001->GetXaxis()->SetTitleOffset(0.9);

  Graph_Graph03001->GetYaxis()->SetLabelFont(42);
  Graph_Graph03001->GetYaxis()->SetLabelOffset(0.004);
  Graph_Graph03001->GetYaxis()->SetLabelSize(0.045);
  Graph_Graph03001->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph03001->GetYaxis()->SetTickLength(0.03);
  Graph_Graph03001->GetYaxis()->SetTitleFont(42);
  Graph_Graph03001->GetYaxis()->SetNdivisions(510, kTRUE);
  Graph_Graph03001->GetYaxis()->SetTitleOffset(1.5);
  //  g->GetHistogram()->Delete();
}

void setpad(TPad *c_1)
{
  // c_1->SetLogx();
  c_1->Draw();
  c_1->cd();
  //   c_1->Range(-0.427002,-0.01611993,4.141898,0.06259125);
  c_1->SetFillColor(0);
  c_1->SetBorderMode(0);
  c_1->SetBorderSize(1);
  c_1->SetTickx(1);
  c_1->SetTicky(1);
  c_1->SetLeftMargin(0.15);
  c_1->SetBottomMargin(0.11);
  c_1->SetTopMargin(0.06);
  c_1->SetRightMargin(0.06);
  c_1->SetFrameBorderMode(0);
  c_1->SetFrameLineWidth(1);
  c_1->SetFrameBorderMode(0);
}

void setlegendstyle(TLegend *leg1)
{
  leg1->SetBorderSize(1);
  leg1->SetLineColor(0);
  leg1->SetLineStyle(1);

  leg1->SetLineWidth(1);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(1001);
}

void scaltepThist(TH1D *h, double scale)
{
  for (int i = 1; i <= h->GetNbinsX(); i++)
  {
    h->SetBinContent(i, h->GetBinContent(i) / h->GetBinWidth(i));
    // h->SetBinError(i,h->GetBinError(i)/h->GetBinWidth(i));
  }
  h->Scale(scale);
}

void setxerrgraphzero(TGraphErrors *g)
{
  for (int i = 0; i < g->GetN(); i++)
  {
    g->SetPointError(i, 0, g->GetErrorY(i));
  }
}

void DivideTGraphErrors(TGraphErrors *g1, TGraphErrors *g2)
{
  int n = g1->GetN();
  for (int i = 0; i < n; i++)
  {
    double x1 = g1->GetPointX(i);
    double x2 = g2->GetPointX(i);
    double y1 = g1->GetPointY(i);
    double y2 = g2->GetPointY(i);

    double x1err = g1->GetErrorX(i);
    double x2err = g2->GetErrorX(i);
    double y1err = g1->GetErrorY(i);
    double y2err = g2->GetErrorY(i);

    double X = x1;
    double Y = y1 / y2;

    double Xerr = x1err;
    double Yerr = Y * TMath::Sqrt(TMath::Power(y1err / y1, 2) +
                                  TMath::Power(y2err / y2, 2));
    g1->SetPoint(i, X, Y);
    g1->SetPointError(i, Xerr, Yerr);
  }
}

const int ndeta = 4;
TString deltaeta[ndeta] = {"", "_eg0.0", "_eg1.0", "_eg2.0"};
TString deltaetaleg[ndeta] = {"", "|#Delta#eta| > 0", "|#Delta#eta| > 1.0",
                              "|#Delta#eta| > 2.0"};
const Char_t deltaetaa[ndeta][10] = {"", "_eg0.0", "_eg1.0", "_eg2.0"};
const Char_t deltaetalegg[ndeta][25] = {"", "|#Delta#eta| > 0", "|#Delta#eta| > 1.0", "|#Delta#eta| > 2.0"};

const int whichetagap = 2; // 2;

Int_t centslice[7] = {0, 10, 20, 30, 50, 70, 100};
const int nflow = 2;
int nflowval[nflow] = {2, 3};

void Plotvncn_NeNe()
{
  TFile *filereadnch = TFile::Open("avgmult_SuperEvent_NeNeWS.root");
  TProfile *pfnch = (TProfile *)filereadnch->Get("avgmult");

  const int nmult = (int)pfnch->GetNbinsX();

  Int_t colors[nflow][ndeta] = {kYellow + 3, kMagenta, kBlue + 1,
                                kRed + 1, kMagenta + 2, kGreen + 2,
                                kAzure + 7, kPink + 5};
  Int_t markers[nflow][ndeta] = {25, 33, 20, 21, 30, 27, 23, 29};
  Double_t size[nflow][ndeta] = {1.8, 1.8, 2.0, 2.0, 1.8, 1.8, 2.0, 2.5};

  TGraphErrors *gcnmultWS[nflow][ndeta];
  TGraphErrors *gvnmultWS[nflow][ndeta];

  TGraphErrors *gdnptWS[nmult][nflow][ndeta];
  TGraphErrors *gvnptWS[nmult][nflow][ndeta];

  TGraphErrors *gcnmultAC[nflow][ndeta];
  TGraphErrors *gvnmultAC[nflow][ndeta];

  TGraphErrors *gvnptratio[nmult][nflow][ndeta];

  TGraphErrors *gdnptAC[nmult][nflow][ndeta];
  TGraphErrors *gvnptAC[nmult][nflow][ndeta];

  TGraphErrors *gcnmultratio[nflow][ndeta];
  TGraphErrors *gvnmultratio[nflow][ndeta];

  TFile *fileWS = TFile::Open("vncn_NeNeWS.root");
  TFile *fileAC = TFile::Open("vncn_NeNeAC.root");

  for (int i = 0; i < nmult; i++)
  {
    TCanvas *c1 = new TCanvas("", "", 750, 800);
    setpad(c1);

    TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
    setpad(pad11);
    pad11->SetBottomMargin(0);
    pad11->Draw();

    TH2D *frame1 = new TH2D("", "", 100, 0.1, 4, 100, -0.03, 0.45);
    sethisto(frame1);
    frame1->SetTitle(";p_{T} (GeV/c);v_{n}{2}(p_{T})");
    frame1->Draw();
    c1->cd(); // Go back to the main canvas before defining pad2
    TPad *pad12 = new TPad("pad12", "pad12", 0, 0.0, 1, 0.3);
    setpad(pad12);
    // pad12->SetLogx();
    // pad12->SetLogy();
    pad12->SetTopMargin(0);
    pad12->SetBottomMargin(0.3);
    pad12->Draw();

    TH2D *hframe12 = new TH2D("", "", 100, 0.1, 4, 100, 0.4, 1.6);
    hframe12->SetTitle(";p_{T} (GeV/c); #frac{Woods-Saxon}{#alpha-cluster}");
    hframe12->GetYaxis()->CenterTitle(true);
    hframe12->SetDirectory(0);
    hframe12->SetStats(0);
    hframe12->GetXaxis()->SetLabelFont(42);
    hframe12->GetXaxis()->SetLabelOffset(0.01);
    hframe12->GetXaxis()->SetLabelSize(0.1);
    hframe12->GetXaxis()->SetTitleSize(0.12);
    hframe12->GetXaxis()->SetTickLength(0.09);
    hframe12->GetXaxis()->SetTitleOffset(1.1);
    hframe12->GetXaxis()->SetTitleFont(42);
    hframe12->GetXaxis()->SetNdivisions(510, kTRUE);
    hframe12->GetYaxis()->SetLabelFont(42);
    hframe12->GetYaxis()->SetLabelOffset(0.01);
    hframe12->GetYaxis()->SetLabelSize(0.1);
    hframe12->GetYaxis()->SetTitleSize(0.1);
    hframe12->GetYaxis()->SetTickLength(0.03);
    hframe12->GetYaxis()->SetTitleFont(42);
    hframe12->GetYaxis()->SetNdivisions(505, kTRUE);
    hframe12->GetYaxis()->SetTitleOffset(0.7);
    hframe12->GetXaxis()->SetMoreLogLabels();
    hframe12->Draw();

    TLine *l1 = new TLine(0.1, 1, 4, 1);
    l1->SetLineStyle(3);
    l1->SetLineColor(kBlack);
    l1->SetLineWidth(2);
    l1->Draw("lsame");

    TCanvas *c2 = new TCanvas("", "", 750, 600);
    setpad(c2);
    TH2D *frame2 = new TH2D("", "", 100, 0.1, 4, 100, -0.005, 0.035);
    sethisto(frame2);
    frame2->SetTitle(";p_{T} (GeV/c);d_{n}{2}(p_{T})");
    frame2->Draw();

    TLegend *leg1 =
        new TLegend(0.212, 0.672, 0.897, 0.884, NULL, "brNDC");
    setlegendstyle(leg1);
    leg1->SetNColumns(2);
    leg1->SetTextSize(0.04);
    leg1->SetHeader("Ne#minusNe, #sqrt{s_{NN}} = 5.36 TeV, IP-Glasma + MUSIC + iSS + UrQMD");
    leg1->AddEntry("", Form("%s", deltaetalegg[whichetagap]), "");
    leg1->AddEntry("", "|#eta| < 0.8", "");
    pad11->cd();
    if (i == 0)
      leg1->Draw();

    TLegend *legcent =
        new TLegend(0.2245989, 0.5563328, 0.3449198, 0.6226068, NULL, "brNDC");
    setlegendstyle(legcent);
    legcent->SetTextSize(0.04);
    legcent->SetHeader(Form("#bf{(%d#minus%d)%%}", centslice[i], centslice[i + 1]));
    legcent->Draw();

    c2->cd();
    legcent->Draw();

    for (int j = 0; j < nflow; j++)
    {
      // for (int k = 1; k < ndeta; k++)
      int k = whichetagap;
      {
        gdnptWS[i][j][k] = (TGraphErrors *)fileWS->Get(
            Form("d%d_2", nflowval[j]) + deltaeta[k] + Form("_pT_mult%d", i));
        gvnptWS[i][j][k] = (TGraphErrors *)fileWS->Get(
            Form("v%d_2", nflowval[j]) + deltaeta[k] + Form("_pT_mult%d", i));

        gdnptAC[i][j][k] = (TGraphErrors *)fileAC->Get(
            Form("d%d_2", nflowval[j]) + deltaeta[k] + Form("_pT_mult%d", i));
        gvnptAC[i][j][k] = (TGraphErrors *)fileAC->Get(
            Form("v%d_2", nflowval[j]) + deltaeta[k] + Form("_pT_mult%d", i));

        c1->cd();
        pad11->cd();
        // setxerrgraphzero(gvnptWS[i][j][k]);
        gvnptWS[i][j][k]->SetLineColor(colors[j][2]);
        gvnptWS[i][j][k]->SetFillStyle(1001);
        gvnptWS[i][j][k]->SetLineWidth(2);
        gvnptWS[i][j][k]->SetLineStyle(1);
        gvnptWS[i][j][k]->SetMarkerStyle(markers[j][2]);
        gvnptWS[i][j][k]->SetMarkerSize(size[j][2]);
        gvnptWS[i][j][k]->SetMarkerColor(colors[j][2]);
        gvnptWS[i][j][k]->SetFillColorAlpha(colors[j][2], 0.8);
        gvnptWS[i][j][k]->Draw("E1,PSAME");

        //  setxerrgraphzero(gvnptAC[i][j][k]);
        gvnptAC[i][j][k]->SetLineColor(colors[j][3]);
        gvnptAC[i][j][k]->SetFillStyle(1001);
        gvnptAC[i][j][k]->SetLineWidth(2);
        gvnptAC[i][j][k]->SetLineStyle(1);
        gvnptAC[i][j][k]->SetMarkerStyle(markers[j][3]);
        gvnptAC[i][j][k]->SetMarkerSize(size[j][3]);
        gvnptAC[i][j][k]->SetMarkerColor(colors[j][3]);
        gvnptAC[i][j][k]->SetFillColorAlpha(colors[j][3], 0.8);
        gvnptAC[i][j][k]->Draw("E1,PSAME");

        c1->cd();
        pad12->cd();
        gvnptratio[i][j][k] = (TGraphErrors *)gvnptWS[i][j][k]->Clone();
        DivideTGraphErrors(gvnptratio[i][j][k], gvnptAC[i][j][k]);
        gvnptratio[i][j][k]->SetLineColor(colors[j][2]);
        gvnptratio[i][j][k]->SetFillStyle(1001);
        gvnptratio[i][j][k]->SetLineWidth(2);
        gvnptratio[i][j][k]->SetLineStyle(1);
        gvnptratio[i][j][k]->SetMarkerStyle(markers[j][2]);
        gvnptratio[i][j][k]->SetMarkerSize(size[j][2]);
        gvnptratio[i][j][k]->SetMarkerColor(colors[j][2]);
        gvnptratio[i][j][k]->SetFillColorAlpha(colors[j][2], 0.8);
        gvnptratio[i][j][k]->Draw("E1,PSAME");

        c2->cd();
        // setxerrgraphzero(gdnptWS[i][j][k]);
        gdnptWS[i][j][k]->SetLineColor(colors[j][2]);
        gdnptWS[i][j][k]->SetFillStyle(1001);
        gdnptWS[i][j][k]->SetLineWidth(2);
        gdnptWS[i][j][k]->SetLineStyle(1);
        gdnptWS[i][j][k]->SetMarkerStyle(markers[j][2]);
        gdnptWS[i][j][k]->SetMarkerSize(size[j][2]);
        gdnptWS[i][j][k]->SetMarkerColor(colors[j][2]);
        gdnptWS[i][j][k]->SetFillColorAlpha(colors[j][2], 0.8);
        gdnptWS[i][j][k]->Draw("E1,PSAME");

        // setxerrgraphzero(gdnptAC[i][j][k]);
        gdnptAC[i][j][k]->SetLineColor(colors[j][3]);
        gdnptAC[i][j][k]->SetFillStyle(1001);
        gdnptAC[i][j][k]->SetLineWidth(2);
        gdnptAC[i][j][k]->SetLineStyle(1);
        gdnptAC[i][j][k]->SetMarkerStyle(markers[j][3]);
        gdnptAC[i][j][k]->SetMarkerSize(size[j][3]);
        gdnptAC[i][j][k]->SetMarkerColor(colors[j][3]);
        gdnptAC[i][j][k]->SetFillColorAlpha(colors[j][3], 0.8);
        gdnptAC[i][j][k]->Draw("E1,PSAME");

        leg1->AddEntry(gvnptWS[i][j][k],
                       Form("n = %d, Woods-Saxon", nflowval[j]), "lp");
        leg1->AddEntry(gvnptAC[i][j][k],
                       Form("n = %d, #alpha-cluster", nflowval[j]), "lp");
      }
    }

    c1->cd();
    c1->SaveAs(Form("vn-mult%d_2_V0MSliceComp%s.pdf", i, deltaetaa[whichetagap]));
    c2->cd();
    leg1->Draw();
    c2->SaveAs(Form("dn-mult%d_2_V0MSliceComp%s.pdf", i, deltaetaa[whichetagap]));
  }

  TCanvas *c1 = new TCanvas("", "", 750, 800);
  setpad(c1);

  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  setpad(pad11);
  pad11->SetBottomMargin(0);
  pad11->Draw();
  // pad11->SetLogx();
  // pad11->SetLogy();

  TH2D *frame11 = new TH2D("", "", 100, 0, 75, 100, -0.005, 0.15);
  sethisto(frame11);
  frame11->SetTitle(";Centrality (%);v_{n}{2} =#sqrt{c_{n}{2}}");
  frame11->Draw();

  TCanvas *c2 = new TCanvas("", "", 750, 800);
  setpad(c2);

  TPad *pad21 = new TPad("pad21", "pad21", 0, 0.3, 1, 1.0);
  setpad(pad21);
  pad21->SetBottomMargin(0);
  pad21->Draw();

  TH2D *frame21 = new TH2D("", "", 100, 0, 75, 100, -0.0005, 0.009);
  sethisto(frame21);
  frame21->SetTitle(";Centrality (%);c_{n}{2}");
  frame21->Draw();

  TLegend *leg1 =
      new TLegend(0.212, 0.672, 0.897, 0.884, NULL, "brNDC");
  setlegendstyle(leg1);
  leg1->SetNColumns(2);
  leg1->SetTextSize(0.04);
  leg1->SetHeader(
      Form("#splitline{Ne#minusNe, #sqrt{s_{NN}} = 5.36 TeV, |#eta| < 0.8, %s}{IP-Glasma +  MUSIC + iSS + UrQMD}", deltaetalegg[whichetagap]));
  pad11->cd();
  leg1->Draw();
  pad21->cd();
  leg1->Draw();

  c1->cd(); // Go back to the main canvas before defining pad2
  TPad *pad12 = new TPad("pad12", "pad12", 0, 0.0, 1, 0.3);
  setpad(pad12);
  // pad12->SetLogx();
  // pad12->SetLogy();
  pad12->SetTopMargin(0);
  pad12->SetBottomMargin(0.3);
  pad12->Draw();

  TH2D *hframe12 = new TH2D("", "", 100, 0, 75, 100, 0.3, 1.7);
  hframe12->SetTitle(";Centrality (%);#frac{Woods-Saxon}{#alpha-cluster}");
  hframe12->GetYaxis()->CenterTitle(true);
  hframe12->SetDirectory(0);
  hframe12->SetStats(0);
  hframe12->GetXaxis()->SetLabelFont(42);
  hframe12->GetXaxis()->SetLabelOffset(0.01);
  hframe12->GetXaxis()->SetLabelSize(0.1);
  hframe12->GetXaxis()->SetTitleSize(0.12);
  hframe12->GetXaxis()->SetTickLength(0.09);
  hframe12->GetXaxis()->SetTitleOffset(1.1);
  hframe12->GetXaxis()->SetTitleFont(42);
  hframe12->GetXaxis()->SetNdivisions(510, kTRUE);
  hframe12->GetYaxis()->SetLabelFont(42);
  hframe12->GetYaxis()->SetLabelOffset(0.01);
  hframe12->GetYaxis()->SetLabelSize(0.1);
  hframe12->GetYaxis()->SetTitleSize(0.1);
  hframe12->GetYaxis()->SetTickLength(0.03);
  hframe12->GetYaxis()->SetTitleFont(42);
  hframe12->GetYaxis()->SetNdivisions(505, kTRUE);
  hframe12->GetYaxis()->SetTitleOffset(0.7);
  hframe12->GetXaxis()->SetMoreLogLabels();
  hframe12->Draw();

  c2->cd(); // Go back to the main canvas before defining pad2
  TPad *pad22 = new TPad("pad22", "pad22", 0, 0.0, 1, 0.3);
  setpad(pad22);
  // pad22->SetLogx();
  // pad22->SetLogy();
  pad22->SetTopMargin(0);
  pad22->SetBottomMargin(0.3);
  pad22->Draw();

  TH2D *hframe22 = new TH2D("", "", 100, 0, 75, 100, 0.05, 2.1);
  hframe22->SetTitle(";Centrality (%);#frac{Woods-Saxon}{#alpha-cluster}");
  hframe22->GetYaxis()->CenterTitle(true);
  hframe22->SetDirectory(0);
  hframe22->SetStats(0);
  hframe22->GetXaxis()->SetLabelFont(42);
  hframe22->GetXaxis()->SetLabelOffset(0.01);
  hframe22->GetXaxis()->SetLabelSize(0.1);
  hframe22->GetXaxis()->SetTitleSize(0.12);
  hframe22->GetXaxis()->SetTickLength(0.09);
  hframe22->GetXaxis()->SetTitleOffset(1.1);
  hframe22->GetXaxis()->SetTitleFont(42);
  hframe22->GetXaxis()->SetNdivisions(510, kTRUE);
  hframe22->GetYaxis()->SetLabelFont(42);
  hframe22->GetYaxis()->SetLabelOffset(0.01);
  hframe22->GetYaxis()->SetLabelSize(0.1);
  hframe22->GetYaxis()->SetTitleSize(0.1);
  hframe22->GetYaxis()->SetTickLength(0.03);
  hframe22->GetYaxis()->SetTitleFont(42);
  hframe22->GetYaxis()->SetNdivisions(505, kTRUE);
  hframe22->GetYaxis()->SetTitleOffset(0.7);
  hframe22->GetXaxis()->SetMoreLogLabels();
  hframe22->Draw();

  TLine *l1 = new TLine(0, 1, 75, 1);
  l1->SetLineStyle(3);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);

  for (int j = 0; j < nflow; j++)
  {
    // for (int k = 1; k < ndeta; k++)
    int k = whichetagap;
    {
      gcnmultWS[j][k] =
          (TGraphErrors *)fileWS->Get(Form("c%d_2", nflowval[j]) + deltaeta[k]);
      gvnmultWS[j][k] =
          (TGraphErrors *)fileWS->Get(Form("v%d_2", nflowval[j]) + deltaeta[k]);

      gcnmultAC[j][k] =
          (TGraphErrors *)fileAC->Get(Form("c%d_2", nflowval[j]) + deltaeta[k]);
      gvnmultAC[j][k] =
          (TGraphErrors *)fileAC->Get(Form("v%d_2", nflowval[j]) + deltaeta[k]);

      pad11->cd();

      setxerrgraphzero(gvnmultWS[j][k]);
      gvnmultWS[j][k]->SetLineColor(colors[j][2]);
      gvnmultWS[j][k]->SetFillStyle(1001);
      gvnmultWS[j][k]->SetLineWidth(2);
      gvnmultWS[j][k]->SetLineStyle(1);
      gvnmultWS[j][k]->SetMarkerStyle(markers[j][2]);
      gvnmultWS[j][k]->SetMarkerSize(size[j][2]);
      gvnmultWS[j][k]->SetMarkerColor(colors[j][2]);
      gvnmultWS[j][k]->SetFillColorAlpha(colors[j][2], 0.8);
      gvnmultWS[j][k]->RemovePoint(5);
      gvnmultWS[j][k]->Draw("E1,CPSAME");

      setxerrgraphzero(gvnmultAC[j][k]);
      gvnmultAC[j][k]->SetLineColor(colors[j][3]);
      gvnmultAC[j][k]->SetFillStyle(1001);
      gvnmultAC[j][k]->SetLineWidth(2);
      gvnmultAC[j][k]->SetLineStyle(1);
      gvnmultAC[j][k]->SetMarkerStyle(markers[j][3]);
      gvnmultAC[j][k]->SetMarkerSize(size[j][3]);
      gvnmultAC[j][k]->SetMarkerColor(colors[j][3]);
      gvnmultAC[j][k]->SetFillColorAlpha(colors[j][3], 0.8);
      gvnmultAC[j][k]->RemovePoint(5);
      gvnmultAC[j][k]->Draw("E1,CPSAME");

      gvnmultratio[j][k] = (TGraphErrors *)gvnmultWS[j][k]->Clone();
      DivideTGraphErrors(gvnmultratio[j][k], gvnmultAC[j][k]);

      pad12->cd();
      l1->Draw("lsame");
      gvnmultratio[j][k]->RemovePoint(5);
      gvnmultratio[j][k]->Draw("E1,CPSAME");

      pad21->cd();

      gcnmultWS[j][k]->SetLineColor(colors[j][2]);
      gcnmultWS[j][k]->SetFillStyle(1001);
      gcnmultWS[j][k]->SetLineWidth(2);
      gcnmultWS[j][k]->SetLineStyle(1);
      gcnmultWS[j][k]->SetMarkerStyle(markers[j][2]);
      gcnmultWS[j][k]->SetMarkerSize(size[j][2]);
      gcnmultWS[j][k]->SetMarkerColor(colors[j][2]);
      gcnmultWS[j][k]->SetFillColorAlpha(colors[j][2], 0.8);
      gcnmultWS[j][k]->Draw("E1,PSAME");

      gcnmultAC[j][k]->SetLineColor(colors[j][3]);
      gcnmultAC[j][k]->SetFillStyle(1001);
      gcnmultAC[j][k]->SetLineWidth(2);
      gcnmultAC[j][k]->SetLineStyle(1);
      gcnmultAC[j][k]->SetMarkerStyle(markers[j][3]);
      gcnmultAC[j][k]->SetMarkerSize(size[j][3]);
      gcnmultAC[j][k]->SetMarkerColor(colors[j][3]);
      gcnmultAC[j][k]->SetFillColorAlpha(colors[j][3], 0.8);
      gcnmultAC[j][k]->Draw("E1,PSAME");

      gcnmultratio[j][k] = (TGraphErrors *)gcnmultWS[j][k]->Clone();
      DivideTGraphErrors(gcnmultratio[j][k], gcnmultAC[j][k]);

      pad22->cd();
      l1->Draw("lsame");
      gcnmultratio[j][k]->Draw("E1,PSAME");

      leg1->AddEntry(gvnmultWS[j][k], Form("n = %d, Woods-Saxon", nflowval[j]),
                     "lp");
      leg1->AddEntry(gvnmultAC[j][k],
                     Form("n = %d, #alpha-cluster", nflowval[j]), "lp");
    }
  }

  c1->cd();
  //		leg1->Draw();
  c1->SaveAs(Form("vn_2_V0MSliceComp%s_NeNe.pdf", deltaetaa[whichetagap]));

  c2->cd();
  // leg1->Draw();
  c2->SaveAs(Form("cn_2_V0MSliceComp%s.pdf", deltaetaa[whichetagap]));

  TCanvas *c3 = new TCanvas("", "", 750, 800);
  setpad(c3);

  TPad *pad31 = new TPad("pad31", "pad31", 0, 0.3, 1, 1.0);
  setpad(pad31);
  pad31->SetBottomMargin(0);
  pad31->Draw();

  TH2D *frame31 = new TH2D("", "", 100, 0, 75, 100, 0.05, 0.8);
  sethisto(frame31);
  frame31->SetTitle(";Centrality (%);v_{3}{2}/v_{2}{2}");
  frame31->Draw();

  c3->cd(); // Go back to the main canvas before defining pad2
  TPad *pad32 = new TPad("pad32", "pad32", 0, 0.0, 1, 0.3);
  setpad(pad32);
  // pad32->SetLogx();
  // pad32->SetLogy();
  pad32->SetTopMargin(0);
  pad32->SetBottomMargin(0.3);
  pad32->Draw();

  TLegend *leg3 =
      new TLegend(0.2927807, 0.7201767, 0.8716578, 0.90243, NULL, "brNDC");
  setlegendstyle(leg3);
  leg3->SetNColumns(2);
  leg3->SetTextSize(0.04);
  leg3->SetHeader(
      Form("#splitline{Ne#minusNe, #sqrt{s_{NN}} = 5.36 TeV, |#eta| < 0.8, %s}{IP-Glasma +  MUSIC + iSS + UrQMD}", deltaetalegg[whichetagap]));

  TH2D *hframe32 = new TH2D("", "", 100, 0, 75, 100, 0.2, 1.7);
  hframe32->SetTitle(";Centrality (%); #frac{Woods-Saxon}{#alpha-cluster}");
  hframe32->GetYaxis()->CenterTitle(true);
  hframe32->SetDirectory(0);
  hframe32->SetStats(0);
  hframe32->GetXaxis()->SetLabelFont(42);
  hframe32->GetXaxis()->SetLabelOffset(0.01);
  hframe32->GetXaxis()->SetLabelSize(0.1);
  hframe32->GetXaxis()->SetTitleSize(0.12);
  hframe32->GetXaxis()->SetTickLength(0.09);
  hframe32->GetXaxis()->SetTitleOffset(1.1);
  hframe32->GetXaxis()->SetTitleFont(42);
  hframe32->GetXaxis()->SetNdivisions(510, kTRUE);
  hframe32->GetYaxis()->SetLabelFont(42);
  hframe32->GetYaxis()->SetLabelOffset(0.01);
  hframe32->GetYaxis()->SetLabelSize(0.1);
  hframe32->GetYaxis()->SetTitleSize(0.1);
  hframe32->GetYaxis()->SetTickLength(0.03);
  hframe32->GetYaxis()->SetTitleFont(42);
  hframe32->GetYaxis()->SetNdivisions(505, kTRUE);
  hframe32->GetYaxis()->SetTitleOffset(0.7);
  hframe32->GetXaxis()->SetMoreLogLabels();
  hframe32->Draw();

  TGraphErrors *gratioWS = (TGraphErrors *)gvnmultWS[1][2]->Clone();
  TGraphErrors *gratioAC = (TGraphErrors *)gvnmultAC[1][2]->Clone();

  gratioWS->SetMarkerColor(kGray + 3);
  gratioWS->SetLineColor(kGray + 3);

  gratioAC->SetMarkerColor(kMagenta);
  gratioAC->SetLineColor(kMagenta);

  DivideTGraphErrors(gratioWS, gvnmultWS[0][2]);
  DivideTGraphErrors(gratioAC, gvnmultAC[0][2]);

  pad31->cd();
  gratioWS->Draw("E1,CPSAME");
  gratioAC->Draw("E1,CPSAME");

  leg3->AddEntry(gratioWS, "Woods-Saxon", "lp");
  leg3->AddEntry(gratioAC, "#alpha-cluster", "lp");

  leg3->Draw();

  TGraphErrors *gratiox = (TGraphErrors *)gratioWS->Clone();
  DivideTGraphErrors(gratiox, gratioAC);

  pad32->cd();
  l1->Draw("lsame");

  gratiox->Draw("E1,CPSAME");

  c3->SaveAs(Form("v2v3ratio%s.pdf", deltaetaa[whichetagap]));
}

Double_t GetPhi(Double_t fpx, Double_t fpy)
{
  const double pi = 3.14159;
  double phi = 0, phi1 = 0;
  // Double_t phi1=0;
  if (fpx == 0 && fpy > 0)
    phi = pi / 2.0;
  if (fpx == 0 && fpy < 0)
    phi = 3 * pi / 2.0;
  if (fpy == 0 && fpx > 0)
    phi = 0;
  if (fpy == 0 && fpx < 0)
    phi = pi;
  phi1 = TMath::ATan(TMath::Abs(fpy / fpx));
  if (fpx > 0 && fpy > 0)
    phi = phi1;
  if (fpx < 0 && fpy > 0)
    phi = pi - phi1;
  if (fpx < 0 && fpy < 0)
    phi = pi + phi1;
  if (fpx > 0 && fpy < 0)
    phi = 2 * pi - phi1;
  return phi;
}

Double_t GetEnergy(Double_t P, Double_t m)
{
  Double_t E;
  E = sqrt(P * P + m * m);
  return E;
}

Double_t GetRapidity(Double_t E, Double_t Pz)
{
  Double_t y;
  if (!(E == Pz))
  {
    y = 0.5 * TMath::Log((E + Pz) / (E - Pz));
  }
  return y;
}

Double_t GetTransMom(Double_t Px, Double_t Py)
{
  Double_t pT;
  pT = sqrt(Px * Px + Py * Py);
  return pT;
}

Double_t GetTotMom(Double_t Px, Double_t Py, Double_t Pz)
{
  Double_t P;
  P = sqrt(Px * Px + Py * Py + Pz * Pz);
  return P;
}

Double_t GetTransMass(Double_t m, Double_t pT)
{
  Double_t mT;
  mT = sqrt(m * m + pT * pT);
  return mT;
}

Double_t GetTheta(Double_t Pz, Double_t P)
{
  Double_t Theta;
  Theta = TMath::ACos(Pz / P);
  return Theta;
}

Double_t GetEta(Double_t theta)
{
  // if(theta == 0) {cout<<"Theta = 0.0; Eta -> inf "<<endl;}
  Double_t Eta;
  Eta = -1.0 * TMath::Log(TMath::Tan(0.5 * theta));
  return Eta;
}

Double_t GetV2(Double_t phi) { return TMath::Cos(2 * phi); }

Double_t GetV3(Double_t phi) { return TMath::Cos(3 * phi); }