/******************************************************************
* Modul name : colormrf.cpp
* Author     : Mihaly Gara (gara@inf.u-szeged.hu) based on the code
*              written by  Csaba Gradwohl (Gradwohl.Csaba@stud.u-szeged.hu)
*              with some  minor contributions from Zoltan Kato
*              (kato@inf.u-szeged.hu).
*
* Copyright  : GNU General Public License www.gnu.org/copyleft/gpl.html
*
* Description:
* Color-based image segmentation using a Markov random field
* segmentation model and four different optimization algorithms:
* Metropolis - Simulated Annealing using Metropolis dynamics
* Gibbs      - Simulated Annealing using a Gibbs sampler
* ICM        - Iterated Conditional Modes, a deterministic suboptimal
*              method (depends on a good initialization).
* MMD        - Modified Metropolis Dynamics, a pseudo-stochastic
*              suboptimal method which is less sensitive to
*              initialization than ICM.
*
* The program GUI is written in wxWidgets hence the code can be
* compiled and ran under Windows as well as under Linux/Unix.
*
* $Id: colormrf.cpp,v 1.1 2009/01/09 20:48:09 kato Exp $
* $Revision: 1.1 $
* $State: Exp $
* $Log: colormrf.cpp,v $
* Revision 1.1  2009/01/09 20:48:09  kato
* Initial revision
*
*
*****************************************************************/
#ifndef lint
static char rcsid_colormrf_cpp[] = "$Id: colormrf.cpp,v 1.1 2009/01/09 20:48:09 kato Exp $";
#endif

/* wxWindows includes
*/
#include <wx/wxprec.h>
#ifndef WX_PRECOMP
#include <wx/wx.h>
#endif
#include <wx/image.h>

#include <wx/colour.h>

#include <math.h>
#include <stdlib.h>
#include <Windows.h>
#include <iostream>
#include <sstream>
#include <cstdlib> 
#include <vector>

/* Random number generators
*/
#include "randomc.h"   // define classes for random number generators

/* Timer classes
*/
#include "CKProcessTimeCounter.h"

#define WINDOW_TITLE "MRF Color Image Segmentation Demo $Revision: 1.1 $"
#define VERSION      "MRF Color Image Segmentation Demo $Revision: 1.1 $" \
  " (Last built " __DATE__" "__TIME__") "
#define COPYRIGHT    "(c) 2006 by Mihaly Gara, Csaba Gradwohl & Zoltan Kato" \
  " (SZTE - Hungary)"

static wxTextCtrl *gaussians;      // output textfield for Gaussian parameters
static CKProcessTimeCounter timer("core"); // CPU timer
static bool timer_valid = FALSE;

/* Program's application class
*/
class MyApp : public wxApp
{
	virtual bool OnInit(); // this is the main entry point
};


/* ImageOperations class: it handles all image operations such as
* loading, saving, etc...
*/
class ImageOperations
{
public:
	ImageOperations(wxWindow *_frame);    // constructor
	wxImage *LoadBmp(wxString bmp_name);	// loads an image from file	
	wxImage *GetOrigImage();
	wxImage *GetLImage();
	wxImage *GetUImage();
	wxImage *GetVImage();
	bool SaveBmp(wxString bmp_name);      // saves out_image to a given file
	bool IsOutput();			// TRUE if  out_image <> NULL
	void SetNoRegions(int n);	      	// sets the number of regions,
	// allocates/frees memory for
	// mean vectors and covariance matrices
	int GetNoRegions() { return no_regions; }
	void SetBeta(double b) { beta = b; }
	void SetT(double x) { t = x; }
	void SetT0(double t) { T0 = t; }
	void SetC(double x) { c = x; }
	void SetAlpha(double x) { alpha = x; }

	//GA
	void SetIterations(double x) { numberOfIterations = (int)x; }
	void SetMProb(double mp) { m_prob = mp; }
	void SetCProb(double cp) { c_prob = cp; }
	void SetChromosomes(double c) { no_chromosomes = (int)c; }

	int GetK() { return K; }
	double GetT() { return T; }
	double GetE() { return E; }
	double GetTimer() { return (timer_valid ? timer.GetElapsedTimeMs() : 0.0); }

	void CalculateMeanAndCovariance(int region);// computes mean and
	// covariance of the given region.
	double CalculateEnergy();                   // computes global energy
	// based on the current
	// lableing in data
	double LocalEnergy(int i, int j, int label);// computes the local
	// energy at site (i,j)
	// assuming "label" has
	// been assigned to it.

	void Metropolis(bool mmd = false);  // executes Metropolis or MMD (if mmd=true)
	void ICM();			    // executes ICM
	void Gibbs();			    // executes Gibbs sampler
	void BeesAlg();           //executes The Bee Algorithm
	void GeneticAlgorithm();



private:
	wxWindow *frame;		    // the main window
	wxImage *in_image, *out_image;    // input & output images. in_image contains
	// the input image in CIE-L*u*v* color space
	wxImage *in_L_image;		    // L* component
	wxImage *in_u_image;		    // u* component
	wxImage *in_v_image;		    // v* component
	int width, height;		    // width and height of the
	// displayed image
	int no_regions;	            // number of regions for Gaussian
	// parameter computation
	int *out_regions;                 // display color of each label (=mean color)
	double beta;                      // strength of second order clique potential
	double t;			    // Stop criteraia threshold: stop
	// if (deltaE < t)
	double T0;		            // Initial temperature (not used by ICM)

	double c;			    // Temperature scheduler's factor:
	// T(n+1)=c*T(n).
	int numberOfIterations;
	double alpha;		            // alpha value for MMD
	double **mean;			    // computed mean values and
	double **variance;		    // variances and 
	double **covariance;		    // covariances for each region
	double **invcov;		    // inverse covariance matrix
	double *denom;                   // denominator for inverse covariance
	double E;			    // current global energy
	double E_old;			    // global energy in the prvious iteration
	double T;			    // current temperature
	int K;			    // current iteration #
	int **classes;		    // this is the labeled image
	double ***in_image_data;	    // Input image (in RGB color space)

	void InitOutImage();
	void SetLuv();		    // Luv settings
	unsigned char *scale(double *luv_vector); // scaling into [0,255]
	double *LuvToRGB(double *luv_pixel);// convert a pixel from CIE-L*u*v* to RGB

	void CreateOutput();	           // creates and draws the output
	// image based on the current labeling
	double Singleton(int i, int inimj, int label); // computes singleton
	// potential at site
	// (i,j) having a label "label"
	double Doubleton(int i, int j, int label); // computes doubleton
	// potential at site
	// (i,j) having a label "label"

	//new variables and methods
	struct Flower {
		int x;
		int y;
		double pollen;
	}; //Flower for our implementation is a given pixcel, so the location, x and y cordination and fitness 
	Flower *scoutWithPollen;

	void placeRandomlyScouts(int n, int start);
	bool scoutExists(Flower flower, int i);
	void evaluateFitness(int n);
	void sortScouts(int n);
	//void selectSiteCenters();

	void ImageOperations::evaluateFitnessPro(int n);
	void siteSearch(int n, int bestSites, int siteSize);

	// genetic algorithm
	struct Chromosome {
		int x;
		int y;
		int val;
		double fitness;
	};

	int **solution; //best solution

	Chromosome **chromosomes; //potential cluster centers
	Chromosome **reordered_chromosomes; //storage for chromosomes during reordering
	Chromosome **mutated_chromosomes; //mutated chromosomes
	std::vector<double> chromosome_energy;
	int no_chromosomes = 20;
	double E_max = 0;
	double E_min = 0;
	double c_prob = .7;
	double m_prob = .05;

	void ImageOperations::makeChromosomes();
	void ImageOperations::evaluate_fitness(int siteSize, bool lastIt);
	void ImageOperations::crossover();
	void ImageOperations::mutation();
	void ImageOperations::formNewPopulation();
	void ImageOperations::printChromosomes();
};


/* MyScrolledWindow class: the window used for diaplaying images */
class MyScrolledWindow : public wxScrolledWindow
{
public:
	MyScrolledWindow(wxWindow* parent, wxWindowID id = -1,
		const wxPoint& pos = wxDefaultPosition,
		const wxSize& size = wxDefaultSize,
		long style = wxHSCROLL | wxVSCROLL,
		const wxString& name = "scrolledWindow") :
		wxScrolledWindow(parent, id, pos, size, style, name)
	{
		bmp = NULL;
	}
	void SetBmp(wxImage *_bmp);     // assigns the image to the window.

protected:
	virtual void OnDraw(wxDC& dc);  // displays the image in the window

private:
	wxImage *bmp;			  // the image to be displayed
	int xDst, yDst;		  // the position of the image within
	// the window (meaningful only when
	// the image is smaller than the window)

	wxMemoryDC memDC;		  // memDC storing the image

	void OnLeftDown(wxMouseEvent& event);	   // Left button event handler
	void OnMouseMotion(wxMouseEvent& event); // mouse motion event handler
	DECLARE_EVENT_TABLE()
};


/* MyFrame class: the main window
*/
class MyFrame : public wxFrame
{
public:
	MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
	~MyFrame();
	/* returns the coordinates of the given training rectangle (or the
	* current one if region==-1)
	*/
	void GetRegion(int &x, int &y, int &w, int &h, int region = -1);

	int GetActRegion() { return act_region; }
	void SetRegs1(int x, int y) {
		regs[act_region * 4] = x;
		regs[act_region * 4 + 1] = y;
	}
	void SetRegs2(int x, int y) {
		regs[act_region * 4 + 2] = x;
		regs[act_region * 4 + 3] = y;
	}
	MyScrolledWindow *GetInputWindow() { return input_window; }
	MyScrolledWindow *GetOutputWindow() { return output_window; }

	bool IsSelected(int region) { // tells whether the region has been selected 
		return (regs[region * 4] || regs[region * 4 + 1] ||
			regs[region * 4 + 2] || regs[region * 4 + 3]);
	}
	bool AllRegionsSelected() { // tells whether all the regions has been selected
		for (int i = 1; i<imageop->GetNoRegions(); ++i)
			if (!IsSelected(i)) return false;
		return true;
	}

private:
	ImageOperations *imageop;
	MyScrolledWindow *input_window, *output_window;    // input & output
	// images' window
	wxChoice *luv_choice;		// scroll-list of CIE-L*u*v* components
	wxButton *load_button, *save_button, *doit_button; // buttons
	wxButton *select_region_button;
	wxChoice *op_choice;		// scroll-list of optimization algorithms
	wxTextCtrl *regions;          // input field for number of classes,
	wxTextCtrl *tbeta, *tt;	// beta, threshold t,
	wxTextCtrl *tT0, *tc;		// initial temperature T0, scheduler factor c,
	wxTextCtrl *talpha;		// and MMD's alpha

	//edited
	wxTextCtrl *number_chromosomes;
	wxTextCtrl *mutation_prob;
	wxTextCtrl *crossover_prob;
	wxTextCtrl *iterations_setter;
	// end of edited

	int act_region;   // the current class
	int *regs;	    // stores the training rectangles for each class.

	/* Event handlers
	*/
	void OnOpen(wxCommandEvent& event);         // Load
	void OnSave(wxCommandEvent& event);         // Save
	void OnDoit(wxCommandEvent& event);         // DoIt
	void OnChoice(wxCommandEvent& event);	      // optimization method selection 
	void OnLuvChoice(wxCommandEvent& event);    // CIE-L*u*v* channel selection
	void OnRegions(wxCommandEvent& event);      // number of classes
	void OnSelectRegion(wxCommandEvent& event); // select training rectangle
	void OnPaint(wxPaintEvent& event);	      // paint handler
	DECLARE_EVENT_TABLE()
};

enum {
	ID_LOAD_BUTTON, ID_SAVE_BUTTON, ID_DOIT_BUTTON, ID_CHOICE, ID_LUV_CHOICE,
	ID_REGIONS, ID_SELECTREGION_BUTTON, ID_BETA, ID_T, ID_T0, ID_C,
	ID_ALPHA, ID_GAUSSIANS, ID_CHROMOSOMES, ID_MUTATIONPROB, ID_CROSSOVERPROB, ID_ITERATIONS
};

/* Event table //runs events on item change?
*/
BEGIN_EVENT_TABLE(MyFrame, wxFrame)
EVT_BUTTON(ID_LOAD_BUTTON, MyFrame::OnOpen)
EVT_BUTTON(ID_SAVE_BUTTON, MyFrame::OnSave)
EVT_BUTTON(ID_DOIT_BUTTON, MyFrame::OnDoit)
EVT_CHOICE(ID_CHOICE, MyFrame::OnChoice)
EVT_CHOICE(ID_LUV_CHOICE, MyFrame::OnLuvChoice)
EVT_PAINT(MyFrame::OnPaint)
EVT_TEXT(ID_REGIONS, MyFrame::OnRegions)
EVT_BUTTON(ID_SELECTREGION_BUTTON, MyFrame::OnSelectRegion)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(MyScrolledWindow, wxScrolledWindow)
EVT_LEFT_DOWN(MyScrolledWindow::OnLeftDown)
EVT_MOTION(MyScrolledWindow::OnMouseMotion)
END_EVENT_TABLE()

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
	MyFrame *frame = new MyFrame(WINDOW_TITLE,
		wxPoint(0, 0), wxSize(1300, 700));
	frame->Show(TRUE);
	frame->Centre(wxBOTH);
	SetTopWindow(frame);
	return TRUE;
}


void MyScrolledWindow::SetBmp(wxImage *_bmp)
{
	xDst = yDst = 0;
	// center if image is smaller than the window
	if (_bmp != NULL)
	{
		if (_bmp->GetWidth() < 300) xDst = (300 - _bmp->GetWidth()) / 2;
		if (_bmp->GetHeight() < 250) yDst = (250 - _bmp->GetHeight()) / 2;
#if wxCHECK_VERSION(2,6,0) // for version 2.6.0 or later
		wxBitmap *mp = new wxBitmap((const wxImage&)*_bmp);
		memDC.SelectObject((wxBitmap&)*mp);
#else    // for version 2.4.x
		memDC.SelectObject(*_bmp);
#endif
	}
	bmp = _bmp;
}


void MyScrolledWindow::OnDraw(wxDC& dc)
{
	if (bmp != NULL)
	{
		// determine which part of the image is visible in the window.
		int x, y;
		GetViewStart(&x, &y);
		x *= 10; y *= 10;	   // must be multiplied by ScrollUnit
		// copy the visible part into the window
		int _xDst, _yDst;
		CalcUnscrolledPosition(xDst, yDst, &_xDst, &_yDst);
		dc.Blit(_xDst, _yDst, 500, 450, &memDC, x, y);

		// draw the training rectangle on the input image
		MyFrame *parent = (MyFrame *)GetParent();
		if (parent->GetInputWindow() == this)
		{
			int x1, y1, w, h;
			parent->GetRegion(x1, y1, w, h);
			if (x1 != 0 || y1 != 0 || w != 0 || h != 0)
			{
				wxPen pen(*wxRED_PEN);
				wxBrush brush(*wxTRANSPARENT_BRUSH);
				dc.SetPen(pen);
				dc.SetBrush(brush);
				dc.DrawRectangle(x1 + xDst, y1 + yDst, w, h);
			}
		}
	}
}


/* Left mouse button event handler
*/
void MyScrolledWindow::OnLeftDown(wxMouseEvent& event)
{
	// comvert window coordinates to image coordinates
	int x, y;
	GetViewStart(&x, &y);
	x *= 10; y *= 10;	 // must be multiplied by ScrollUnit

	MyFrame *frame = (MyFrame *)GetParent();
	if (frame->GetActRegion() != -1)	// in this case regs != NULL
	{
		if (event.m_x >= xDst && event.m_x < xDst + bmp->GetWidth() &&
			event.m_y >= yDst && event.m_y < yDst + bmp->GetHeight())
			frame->SetRegs1(event.m_x + x - xDst, event.m_y + y - yDst); // scroll added
	}
}


/* Mouse motion event handler
*/
void MyScrolledWindow::OnMouseMotion(wxMouseEvent& event)
{
	if (event.LeftIsDown())
	{
		// comvert window coordinates to image coordinates
		int x, y;
		GetViewStart(&x, &y);
		x *= 10; y *= 10;		// must be multiplied by ScrollUnit
		MyFrame *frame = (MyFrame *)GetParent();
		if (frame->GetActRegion() != -1)		// in this case regs != NULL
		{
			if (event.m_x >= xDst && event.m_x < xDst + bmp->GetWidth() &&
				event.m_y >= yDst && event.m_y < yDst + bmp->GetHeight())
				frame->SetRegs2(event.m_x + x - xDst, event.m_y + y - yDst); // scroll added
			Refresh();
		}
	}
}


/*********************************************************************
/* Functions of MyFrame class
/********************************************************************/
void MyFrame::OnPaint(wxPaintEvent& event)
{
	wxPaintDC pDC(this);

	wxString str;
	str.Printf("Number of classes:");
	pDC.DrawText(str, 20, 525);

	str.Printf("ß = ");
	pDC.DrawText(str, 43, 560);
	str.Printf("t = ");
	pDC.DrawText(str, 181, 560);
	str.Printf("Class parameters:");
	pDC.DrawText(str, 20, 665);
	if (op_choice->GetStringSelection() != "ICM")
	{
		str.Printf("T0 = ");
		pDC.DrawText(str, 35, 595);
		str.Printf("c = ");
		pDC.DrawText(str, 181, 595);
	}
	str.Printf(VERSION);
	pDC.DrawText(str, 20, 790);
	str.Printf(COPYRIGHT);
	pDC.DrawText(str, 20, 805);
	//  str.Printf(ADD_COLOR); //??
	//  pDC.DrawText(str, 20, 620);

	//edited 

	if (op_choice->GetStringSelection() == "GeneticAlgorithm") {

		str.Printf("# of Iterations:");
		pDC.DrawText(str, 340, 525);

		str.Printf("Mutation Probability:");
		pDC.DrawText(str, 300, 560);

		str.Printf("Crossover Probability:");
		pDC.DrawText(str, 290, 595);

		str.Printf("# of Chromosomes:");
		pDC.DrawText(str, 310, 630);
	}


	//end edited

	if (op_choice->GetStringSelection() == "MMD")
	{
		str.Printf("alpha = ");
		pDC.DrawText(str, 15, 630);
	}

	str.Printf("iteration = ");
	pDC.DrawText(str, 675, 560);
	str.Printf("global energy = ");
	pDC.DrawText(str, 675, 595);
	str.Printf("T = ");
	pDC.DrawText(str, 675, 630);
	str.Printf("CPU time = ");
	pDC.DrawText(str, 675, 665);
	str.Printf("Height");
	pDC.DrawText(str, 675, 700);
	str.Printf("Width");
	pDC.DrawText(str, 675, 735);
	pDC.DrawText(wxString() << imageop->GetK(), 805, 560);
	pDC.DrawText(wxString() << imageop->GetE(), 805, 595);
	pDC.DrawText(wxString() << imageop->GetT(), 805, 630);
	pDC.DrawText(wxString() << imageop->GetTimer() << " ms", 805, 665);
	event.Skip();
}


MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
	imageop = new ImageOperations(this);
	input_window = new MyScrolledWindow(this, -1, wxPoint(15, 15),
		wxSize(500, 450), wxHSCROLL | wxVSCROLL);
	input_window->SetBackgroundColour(wxColour(255, 255, 255));
	output_window = new MyScrolledWindow(this, -1, wxPoint(675, 15),
		wxSize(500, 450), wxHSCROLL | wxVSCROLL);
	output_window->SetBackgroundColour(wxColour(255, 255, 255));

	wxString luv_label[4] = { "Original", "L", "u", "v" };
	luv_choice = new wxChoice(this, ID_LUV_CHOICE, wxPoint(15, 480), wxDefaultSize,
		4, luv_label);
	luv_choice->SetStringSelection("Original");
	luv_choice->Disable();

	load_button = new wxButton(this, ID_LOAD_BUTTON, "Load", wxPoint(125, 480));
	save_button = new wxButton(this, ID_SAVE_BUTTON, "Save", wxPoint(675, 480));
	doit_button = new wxButton(this, ID_DOIT_BUTTON, "Do it >>",
		wxPoint(558, 150));
	doit_button->Disable();
	select_region_button = new wxButton(this, ID_SELECTREGION_BUTTON,
		"Select classes", wxPoint(203, 521));
	select_region_button->Disable();
	wxString choices[6] = { "Metropolis", "Gibbs sampler", "ICM", "MMD", "BeesAlg", "GeneticAlgorithm" };
	op_choice = new wxChoice(this, ID_CHOICE, wxPoint(546, 80), wxDefaultSize,
		6, choices);
	op_choice->SetStringSelection("Metropolis");

	regions = new wxTextCtrl(this, ID_REGIONS, "", wxPoint(152, 521),
		wxSize(27, 20),
		wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	regions->SetMaxLength(3);
	regions->Disable();
	tbeta = new wxTextCtrl(this, ID_BETA, "", wxPoint(67, 556),
		wxSize(60, 20),
		wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	tbeta->SetMaxLength(8);
	//  tbeta->SetValue("2.5");
	*tbeta << 2.5;
	tt = new wxTextCtrl(this, ID_T, "", wxPoint(203, 556),
		wxSize(60, 20), wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	tt->SetMaxLength(8);
	// tt->SetValue("0.05");
	*tt << 0.05;
	tT0 = new wxTextCtrl(this, ID_T0, "", wxPoint(67, 591), wxSize(60, 20),
		wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	tT0->SetMaxLength(8);
	// tT0->SetValue("4.0");
	*tT0 << 4.0;
	tc = new wxTextCtrl(this, ID_C, "", wxPoint(203, 591), wxSize(60, 20),
		wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	tc->SetMaxLength(8);
	// tc->SetValue("0.98");
	*tc << 0.98;
	talpha = new wxTextCtrl(this, ID_ALPHA, "", wxPoint(67, 626),
		wxSize(60, 20), wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	talpha->SetMaxLength(8);
	// talpha->SetValue("0.1");
	*talpha << 0.1;
	talpha->Hide();

	//edited

	mutation_prob = new wxTextCtrl(this, ID_MUTATIONPROB, "", wxPoint(440, 556),
		wxSize(60, 20), wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	mutation_prob->SetMaxLength(8);
	*mutation_prob << 0.1;
	mutation_prob->Hide();

	crossover_prob = new wxTextCtrl(this, ID_CROSSOVERPROB, "", wxPoint(440, 591),
		wxSize(60, 20), wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	crossover_prob->SetMaxLength(8);
	*crossover_prob << 0.1;
	crossover_prob->Hide();

	number_chromosomes = new wxTextCtrl(this, ID_CHROMOSOMES, "", wxPoint(440, 626),
		wxSize(60, 20), wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	number_chromosomes->SetMaxLength(8);
	*number_chromosomes << 100;
	number_chromosomes->Hide();

	iterations_setter = new wxTextCtrl(this, ID_ITERATIONS, "", wxPoint(440, 521),
		wxSize(60, 20), wxTE_PROCESS_ENTER | wxTE_RIGHT,
		*(new wxTextValidator(wxFILTER_NUMERIC)));
	iterations_setter->SetMaxLength(8);
	*iterations_setter << 50;
	iterations_setter->Hide();

	//end edited

	gaussians = new wxTextCtrl(this, ID_GAUSSIANS, "", wxPoint(20, 680),
		wxSize(500, 100),
		wxTE_MULTILINE | wxTE_DONTWRAP | wxTE_READONLY,
		wxDefaultValidator);
	gaussians->SetValue("# Mean (L, u, v)\t\tVariance (L, u, v)\t\tCovariance (L-u, L-v, u-v)\n");

	regs = NULL;
	act_region = -1;

}


MyFrame::~MyFrame()
{
	delete imageop;
}


void MyFrame::OnOpen(wxCommandEvent& event)
{
	wxString image_name;
	wxFileDialog* fdialog = new wxFileDialog(this, "Open file", "", "",
		"BMP files (*.bmp)|*.bmp",
		wxOPEN | wxCHANGE_DIR);

	if (fdialog->ShowModal() == wxID_OK)
	{
		image_name = fdialog->GetPath();
		wxImage *bmp;
		if (bmp = imageop->LoadBmp(image_name)) // image succesfully loaded
		{
			input_window->SetScrollbars(10, 10, (bmp->GetWidth()) / 10,
				(bmp->GetHeight()) / 10);
			input_window->SetBmp(bmp);
			output_window->SetBmp(NULL);
			output_window->SetScrollbars(10, 10, 0, 0);
			input_window->Refresh();
			output_window->Refresh();
			// enable input fields
			regions->Enable();
			luv_choice->Enable();
			luv_choice->SetStringSelection("Original");
			select_region_button->SetLabel("Select classes");  // reset button
			// label
			if (regs != NULL)
			{
				delete[] regs;  // remove all rectangle selections
				regs = NULL;
				act_region = -1;
				imageop->SetNoRegions(-1);
			}
			doit_button->Disable();
			Refresh();
		}
	}
}


void MyFrame::OnSave(wxCommandEvent& event)
{
	if (imageop->IsOutput()) // if there is anything to save
	{
		wxString image_name;
		wxFileDialog* fdialog = new wxFileDialog(this, "Save file as", "", "",
			"BMP files (*.bmp)|*.bmp",
			wxSAVE | wxCHANGE_DIR |
			wxOVERWRITE_PROMPT);

		if (fdialog->ShowModal() == wxID_OK)
		{
			image_name = fdialog->GetPath();
			if (!imageop->SaveBmp(image_name)) // saving failed
				wxLogError("Can't save image!", "ERROR");
			//	    wxMessageBox("Can't save image!", "ERROR");
		}
	}
}


void MyFrame::OnDoit(wxCommandEvent& event)
{

	wxString beta, t, T0, c, alpha;
	wxString iterations, crossoverPb, mutationPb, chrmsmeCount;

	if ((beta = tbeta->GetValue()).Length() == 0)
	{
		wxLogWarning("ß value missing!", "Warning!");
		//     wxMessageBox("ß value missing", "Error");
		return;
	}
	else	// TODO: check value!
		imageop->SetBeta(atof(beta));
	if ((t = tt->GetValue()).Length() == 0)
	{
		wxLogWarning("t value missing!", "Warning!");
		//      wxMessageBox("t value missing", "Error");
		return;
	}
	else	// TODO: check value!
		imageop->SetT(atof(t));

	if (op_choice->GetStringSelection() != "ICM")
	{
		if ((T0 = tT0->GetValue()).Length() == 0)
		{
			wxLogWarning("T0 value missing!", "Warning!");
			//	  wxMessageBox("T0 value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetT0(atof(T0));
		if ((c = tc->GetValue()).Length() == 0)
		{
			wxLogWarning("c value missing!", "Warning!");
			//	  wxMessageBox("c value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetC(atof(c));
	}
	if (op_choice->GetStringSelection() == "MMD")
	{
		if ((alpha = talpha->GetValue()).Length() == 0)
		{
			wxLogWarning("alpha value missing!", "Warning!");
			//	  wxMessageBox("alpha value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetAlpha(atof(alpha));
	}
	if (op_choice->GetStringSelection() == "GeneticAlgorithm")
	{
		//check iterations
		if ((iterations = iterations_setter->GetValue()).Length() == 0)
		{
			wxLogWarning("iteration value missing!", "Warning!");
			//	  wxMessageBox("iteration value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetIterations(atof(iterations));

		//check crossover prob
		if ((crossoverPb = crossover_prob->GetValue()).Length() == 0)
		{
			wxLogWarning("crossover value missing!", "Warning!");
			//	  wxMessageBox("crossover value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetCProb(atof(crossoverPb));

		//check mutation prob
		if ((mutationPb = mutation_prob->GetValue()).Length() == 0)
		{
			wxLogWarning("mutation value missing!", "Warning!");
			//	  wxMessageBox("mutation value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetMProb(atof(mutationPb));

		//check no chromosomes
		if ((chrmsmeCount = number_chromosomes->GetValue()).Length() == 0)
		{
			wxLogWarning("# of chromosomes value missing!", "Warning!");
			//	  wxMessageBox("# of chromosomes value missing", "Error");
			return;
		}
		else	// TODO: check value!
			imageop->SetChromosomes(atof(chrmsmeCount));
	}

	timer_valid = FALSE; // timer's value is invalid. Used by GetTimer()
	Refresh();
	timer.Reset();       // reset timer
	timer.Start();       // start timer
	if (op_choice->GetStringSelection() == "Metropolis")
	{
		imageop->Metropolis();
	}
	else if (op_choice->GetStringSelection() == "MMD")
	{
		imageop->Metropolis(true);
	}
	else if (op_choice->GetStringSelection() == "ICM")
	{
		imageop->ICM();
	}
	else if (op_choice->GetStringSelection() == "Gibbs sampler")
	{
		imageop->Gibbs();
	}
	else if (op_choice->GetStringSelection() == "BeesAlg")
	{
		imageop->BeesAlg();
	}
	else if (op_choice->GetStringSelection() == "GeneticAlgorithm")
	{
		imageop->GeneticAlgorithm();
	}

	timer.Stop();       // stop timer
	timer_valid = TRUE; // timer's value is valid. Used by GetTimer()
	Refresh();
}


/* selection list handler
*/
void MyFrame::OnChoice(wxCommandEvent& event)
{
	if (op_choice->GetStringSelection() == "ICM")
	{
		tT0->Hide();
		tc->Hide();
	}
	else
	{
		tT0->Show();
		tc->Show();
	}
	if (op_choice->GetStringSelection() == "MMD")
	{
		talpha->Show();
	}
	else
	{
		talpha->Hide();
	}

	if (op_choice->GetStringSelection() == "GeneticAlgorithm")
	{
		iterations_setter->Show();
		mutation_prob->Show();
		crossover_prob->Show();
		number_chromosomes->Show();
	}
	else
	{
		iterations_setter->Hide();
		mutation_prob->Hide();
		crossover_prob->Hide();
		number_chromosomes->Hide();
	}
	Refresh();
}

void MyFrame::OnLuvChoice(wxCommandEvent& event)
{
	if (luv_choice->GetStringSelection() == "Original")
	{
		input_window->SetBmp(imageop->GetOrigImage());
	}
	if (luv_choice->GetStringSelection() == "L")
	{
		input_window->SetBmp(imageop->GetLImage());
	}
	if (luv_choice->GetStringSelection() == "u")
	{
		input_window->SetBmp(imageop->GetUImage());
	}
	if (luv_choice->GetStringSelection() == "v")
	{
		input_window->SetBmp(imageop->GetVImage());
	}

	input_window->Refresh();
}

/* called whenever number of regions has been changed.
*/
void MyFrame::OnRegions(wxCommandEvent& event)
{
	select_region_button->SetLabel("Select classes"); // reset button label
	{
		delete[] regs;	// remove all rectangle selections
		regs = NULL;
		act_region = -1;
		imageop->SetNoRegions(-1);
		gaussians->SetValue("# Mean (L, u, v)\t\tVariance (L, u, v)\t\tCovariance (L-u, L-v, u-v)\n");
		// clear parameter textfield
	}
	doit_button->Disable();
	wxString str = regions->GetValue(); // get entered value
	// TODO: check wheter the value is a positive integer!
	if (str.Length() != 0) select_region_button->Enable();
	else select_region_button->Disable();

}


/* on  "Select classes" button pushed
*/
void MyFrame::OnSelectRegion(wxCommandEvent& event)
{
	static int no_regions = 0;
	if (select_region_button->GetLabel() == "Select classes")
	{
		act_region = 0;
		no_regions = atoi(regions->GetValue());     // get number of regions
		imageop->SetNoRegions(no_regions);
		regs = new int[no_regions * 4];                   // aloccate memory
		for (int i = 0; i<no_regions * 4; ++i) regs[i] = 0; // init with 0
		select_region_button->SetLabel(act_region == no_regions - 1 ?
			"Finish" : "Next class");  // change label
	}
	else	// select next training rectangle
	{
		if (act_region == no_regions - 2)
			select_region_button->SetLabel("Finish");  // change label
		if (IsSelected(act_region))
			imageop->CalculateMeanAndCovariance(act_region);
		if (++act_region >= imageop->GetNoRegions())   // no more regions
			act_region = -1;                             // no act_region
		if (AllRegionsSelected()) {
			doit_button->Enable();            // Enable DoIt button
			select_region_button->Disable();  // Disable "Next class" button
		}

	}
	input_window->Refresh();
}


/* return the coordinates of the given training rectangle
*/
void MyFrame::GetRegion(int &x, int &y, int &w, int &h, int region)
{
	x = y = w = h = 0;
	if (region == -1) // -1 ==> return current rectangle
		region = act_region;
	if (region != -1) // otherwise compute coordinates
	{
		int x1 = regs[region * 4];
		int y1 = regs[region * 4 + 1];
		int x2 = regs[region * 4 + 2];
		int y2 = regs[region * 4 + 3];
		x = x1<x2 ? x1 : x2;		// left X coordinate
		y = y1<y2 ? y1 : y2;		// lower Y coordinate
		w = abs(x1 - x2);			// width
		h = abs(y1 - y2);			// height
	}
}


/*********************************************************************
/* Functions of ImageOperations class
/********************************************************************/

//state of application upon initialization
ImageOperations::ImageOperations(wxWindow *_frame)
{
	frame = _frame;
	frame->SetBackgroundColour(wxColour(204, 236, 255));
	in_image = out_image = NULL;

	no_regions = -1; // -1 ==> num. of regions has not been specified yet!
	beta = -1;
	T0 = -1;
	c = -1;
	K = 0;
	E = 0;
	T = 0;
	mean = variance = NULL;
	covariance = invcov = NULL;
	denom = NULL;
	alpha = 0.1;

	no_chromosomes = 20;
	c_prob = .7;
	m_prob = .05;

}


wxImage *ImageOperations::LoadBmp(wxString bmp_name)
{
	wxImage *img = new wxImage(bmp_name);
	if (img->Ok()) // set new values						
	{
		in_image = img;
		height = in_image->GetHeight();
		width = in_image->GetWidth();
		SetLuv();
		out_image = NULL;
	}
	return in_image;
}

wxImage *ImageOperations::GetOrigImage()
{
	return in_image;
}

wxImage *ImageOperations::GetLImage()
{
	return in_L_image;
}

wxImage *ImageOperations::GetUImage()
{
	return in_u_image;
}

wxImage *ImageOperations::GetVImage()
{
	return in_v_image;
}

bool ImageOperations::SaveBmp(wxString bmp_name)
{
	return out_image->SaveFile(bmp_name, wxBITMAP_TYPE_BMP);
}


bool ImageOperations::IsOutput()
{
	if (out_image == NULL) return false;
	else return true;
}


void ImageOperations::SetNoRegions(int n)
{
	int j;
	no_regions = n;
	if (n == -1)
	{
		delete mean;
		delete variance;
		delete covariance;
		delete invcov;
		delete denom;
		mean = variance = NULL;
		covariance = invcov = NULL;
		denom = NULL;
	}
	else
	{
		//2d arrays
		mean = new double*[3];
		variance = new double*[3];
		covariance = new double*[3];
		invcov = new double*[6];
		denom = new double[n];
		for (j = 0; j<3; j++)
		{
			mean[j] = new double[n];
			variance[j] = new double[n];
			covariance[j] = new double[n];
			invcov[j] = new double[n];
		}
		for (j = 3; j<6; j++)
			invcov[j] = new double[n];
		for (int i = 0; i<n; ++i)
		{
			for (j = 0; j<3; j++)
				mean[j][i] = variance[j][i] = covariance[j][i] = invcov[j][i] = -1;
			for (j = 3; j<6; j++)
				invcov[j][i] = -1;
		}
	}
}


/* Compute mean vector and covariance matrix for a given region
*/
void ImageOperations::CalculateMeanAndCovariance(int region)
{
	if (in_image != NULL)
	{
		int x, y, w, h;
		int i, j, k;
		double sum, sum2, sum3;
		((MyFrame *)frame)->GetRegion(x, y, w, h, region);

		for (k = 0; k<3; k++)
		{
			sum = sum2 = 0;
			for (i = y; i<y + h; ++i)
				for (j = x; j<x + w; ++j)
				{
					sum += in_image_data[i][j][k];
					sum2 += in_image_data[i][j][k] * in_image_data[i][j][k];
				}

			mean[k][region] = sum / (w*h);
			variance[k][region] = (sum2 - (sum*sum) / (w*h)) / (w*h - 1);
		}

		// compute covariances
		sum = sum2 = sum3 = 0;
		for (i = y; i<y + h; ++i)
			for (j = x; j<x + w; ++j)
			{		// L-u covariance
				sum += (in_image_data[i][j][0] - mean[0][region])*(in_image_data[i][j][1] - mean[1][region]);
				// L-v covariance
				sum2 += (in_image_data[i][j][0] - mean[0][region])*(in_image_data[i][j][2] - mean[2][region]);
				// u-v covariance
				sum3 += (in_image_data[i][j][1] - mean[1][region])*(in_image_data[i][j][2] - mean[2][region]);
			}

		covariance[0][region] = sum / (w*h);   // L-u covariance
		covariance[1][region] = sum2 / (w*h);  // L-v covariance
		covariance[2][region] = sum3 / (w*h);  // u-v covariance
		// Compute elements of inverse covariance matrix
		// element (1,1)
		invcov[0][region] = variance[2][region] * variance[1][region] - covariance[2][region] * covariance[2][region];

		// elements (1,2) and (2,1)
		invcov[1][region] = covariance[1][region] * covariance[2][region] - variance[2][region] * covariance[0][region];

		// elements (1,3) and (3,1)
		invcov[2][region] = covariance[0][region] * covariance[2][region] - variance[1][region] * covariance[1][region];

		// element (2,2)
		invcov[3][region] = variance[2][region] * variance[0][region] - covariance[1][region] * covariance[1][region];

		// elements (2,3) and (3,2)
		invcov[4][region] = covariance[0][region] * covariance[1][region] - variance[0][region] * covariance[2][region];

		// element (3,3)
		invcov[5][region] = variance[1][region] * variance[0][region] - covariance[0][region] * covariance[0][region];

		// denominator for computing elements of 
		// inverse covariance matrix
		denom[region] = variance[0][region] * variance[1][region] * variance[2][region] -
			variance[2][region] * covariance[0][region] * covariance[0][region] -
			variance[1][region] * covariance[1][region] * covariance[1][region] -
			variance[0][region] * covariance[2][region] * covariance[2][region] +
			covariance[0][region] * covariance[1][region] * covariance[2][region] * 2;

		if (denom[region] == 0)
			denom[region] = 1e-10;
		for (k = 0; k<3; k++)
		{
			if (covariance[k][region] == 0)
				covariance[k][region] = 1e-10;
			if (variance[k][region] == 0)
				variance[k][region] = 1e-10;
		}

		// print parameters in gaussians textfield
		*gaussians << region + 1 << " (" << mean[0][region] << ", " << mean[1][region] << ", " <<
			mean[2][region] << ")\t(" << variance[0][region] << ", " << variance[1][region] << ", "
			<< variance[2][region] << ")\t(" << covariance[0][region] << ", " <<
			covariance[1][region] << ", " << covariance[2][region] << ")\n";
	}
}

double ImageOperations::Singleton(int i, int j, int label)
{
	double det;    // determinant of covariance matrix
	double gauss;  // exponential part of Gaussians

	det = variance[0][label] * variance[1][label] * variance[2][label] +
		2 * covariance[0][label] * covariance[1][label] * covariance[0][label] -
		covariance[0][label] * covariance[0][label] * variance[2][label] -
		covariance[1][label] * covariance[1][label] * variance[1][label] -
		covariance[2][label] * covariance[2][label] * variance[0][label];

	gauss = ((in_image_data[i][j][0] - mean[0][label]) * invcov[0][label] +
		(in_image_data[i][j][1] - mean[1][label]) * invcov[1][label] +
		(in_image_data[i][j][2] - mean[2][label]) * invcov[2][label]) * (in_image_data[i][j][0] - mean[0][label]) +
		((in_image_data[i][j][0] - mean[0][label]) * invcov[1][label] +
		(in_image_data[i][j][1] - mean[1][label]) * invcov[3][label] +
		(in_image_data[i][j][2] - mean[2][label]) * invcov[4][label]) * (in_image_data[i][j][1] - mean[1][label]) +
		((in_image_data[i][j][0] - mean[0][label]) * invcov[2][label] +
		(in_image_data[i][j][1] - mean[1][label]) * invcov[4][label] +
		(in_image_data[i][j][2] - mean[2][label]) * invcov[5][label]) * (in_image_data[i][j][2] - mean[2][label]);

	if (det == 0)
		det = 1e-10;
	else if (det<0)
	{
		det = -det;
		//	  return - log(sqrt(2.0*3.141592653589793*det)) + 0.5 * (double)gauss / (double)denom[label];
	}
	return log(sqrt(2.0*3.141592653589793*det)) + 0.5 * (double)gauss / (double)denom[label];
}


double ImageOperations::Doubleton(int i, int j, int label)
{
	double energy = 0.0;

	if (i != height - 1) // south
	{
		if (label == classes[i + 1][j]) energy -= beta;
		else energy += beta;
	}
	if (j != width - 1) // east
	{
		if (label == classes[i][j + 1]) energy -= beta;
		else energy += beta;
	}
	if (i != 0) // nord
	{
		if (label == classes[i - 1][j]) energy -= beta;
		else energy += beta;
	}
	if (j != 0) // west
	{
		if (label == classes[i][j - 1]) energy -= beta;
		else energy += beta;
	}
	return energy;
}


/* compute global energy
*/
double ImageOperations::CalculateEnergy()
{
	double singletons = 0.0;
	double doubletons = 0.0;
	int i, j, k;
	for (i = 0; i<height; ++i)
		for (j = 0; j<width; ++j)
		{
			k = classes[i][j];
			// singleton
			singletons += Singleton(i, j, k);
			// doubleton
			doubletons += Doubleton(i, j, k); // Note: here each doubleton is
			// counted twice ==> divide by
			// 2 at the end!
		}
	return singletons + doubletons / 2;
}



double ImageOperations::LocalEnergy(int i, int j, int label)
{
	return Singleton(i, j, label) + Doubleton(i, j, label);
}


/* Initialize segmentation
*/
void ImageOperations::InitOutImage()
{
	int i, j, k, r;
	double *temp_data;
	double e, e2;	 // store local energy

	classes = new int*[height]; // allocate memory for classes
	for (i = 0; i<height; ++i)
		classes[i] = new int[width];
	/* initialize using Maximum Likelihood (~ max. of singleton energy)
	*/
	for (i = 0; i<height; ++i)
		for (j = 0; j<width; ++j)
		{
			//local energy
			e = Singleton(i, j, 0);
			classes[i][j] = 0;
			for (r = 1; r<no_regions; ++r)
				//energy minimization
				if ((e2 = Singleton(i, j, r)) < e)
				{
					e = e2;
					classes[i][j] = r;
				}
		}
	out_regions = new int[no_regions * 3];
	temp_data = new double[3];
	for (r = 0; r<no_regions; r++)
	{
		temp_data[0] = mean[0][r];
		temp_data[1] = mean[1][r];
		temp_data[2] = mean[2][r];
		temp_data = LuvToRGB(temp_data);
		for (k = 0; k<3; k++)
		{
			out_regions[r * 3 + k] = (int)temp_data[k];
		}
	}
}

/* Compute CIE-L*u*v* values and
* set in_image, in_L_image, in_u_image, in_v_image
*/
void ImageOperations::SetLuv()
{
	int i, j;
	double *luv_data;
	unsigned char *scaled_luv_data;
	unsigned char *l_data;
	unsigned char *u_data;
	unsigned char *v_data;
	double *xyz_data;
	double u0, v0;

	unsigned char *in_data = in_image->GetData();

	luv_data = (double *)malloc(width*height * 3 * sizeof(double));
	l_data = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));
	u_data = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));
	v_data = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));
	scaled_luv_data = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));
	xyz_data = (double *)malloc(width*height * 3 * sizeof(double));

	// Compute u0, v0 (corresponding to white color)
	u0 = 4 * 242.36628 / (242.36628 + 15 * 254.999745 + 3 * 277.63227);
	v0 = 9 * 254.999754 / (242.36628 + 15 * 254.999745 + 3 * 277.63227);

	// Convert into CIE-XYZ color space
	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			// X component
			xyz_data[(i*width * 3) + j * 3] = (in_data[i*width * 3 + j * 3] * 0.412453 +
				in_data[i*width * 3 + j * 3 + 1] * 0.35758 +
				in_data[i*width * 3 + j * 3 + 2] * 0.180423);
			// Y component
			xyz_data[(i*width * 3) + j * 3 + 1] = (in_data[i*width * 3 + j * 3] * 0.212671 +
				in_data[i*width * 3 + j * 3 + 1] * 0.715160 +
				in_data[i*width * 3 + j * 3 + 2] * 0.072169);
			// Z component
			xyz_data[(i*width * 3) + j * 3 + 2] = (in_data[i*width * 3 + j * 3] * 0.019334 +
				in_data[i*width * 3 + j * 3 + 1] * 0.119193 +
				in_data[i*width * 3 + j * 3 + 2] * 0.950227);
		}

	// Convert into CIE-L*u*v* color space
	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			// L component
			if ((xyz_data[(i*width * 3) + j * 3 + 1] / 254.999745) > 0.008856)
				luv_data[(i*width * 3) + j * 3] = 116 * pow((xyz_data[(i*width * 3) + j * 3 + 1] / 254.999745), (1.0 / 3.0)) - 16;
			else
				luv_data[(i*width * 3) + j * 3] = 903.3*(xyz_data[(i*width * 3) + j * 3 + 1] / 254.999745);

			// u component
			if ((xyz_data[(i*width * 3) + j * 3] + 15 * xyz_data[(i*width * 3) + j * 3 + 1] + 3 * xyz_data[(i*width * 3) + j * 3 + 2]) == 0)
				luv_data[(i*width * 3) + j * 3 + 1] = 13 * luv_data[(i*width * 3) + j * 3] * (-u0);
			else
				luv_data[(i*width * 3) + j * 3 + 1] = 13 * luv_data[(i*width * 3) + j * 3] * ((4 * xyz_data[(i*width * 3) + j * 3] /
				(xyz_data[(i*width * 3) + j * 3] + 15 * xyz_data[(i*width * 3) + j * 3 + 1] + 3 * xyz_data[(i*width * 3) + j * 3 + 2])) - u0);

			// v component
			if ((xyz_data[(i*width * 3) + j * 3] + 15 * xyz_data[(i*width * 3) + j * 3 + 1] + 3 * xyz_data[(i*width * 3) + j * 3 + 2]) == 0)
				luv_data[(i*width * 3) + j * 3 + 2] = 13 * luv_data[(i*width * 3) + j * 3] * (-v0);
			else
				luv_data[(i*width * 3) + j * 3 + 2] = 13 * luv_data[(i*width * 3) + j * 3] * ((9 * xyz_data[(i*width * 3) + j * 3 + 1] /
				(xyz_data[(i*width * 3) + j * 3] + 15 * xyz_data[(i*width * 3) + j * 3 + 1] + 3 * xyz_data[(i*width * 3) + j * 3 + 2])) - v0);
		}

	in_image_data = new double**[height]; // allocate memory for in_image_data
	for (i = 0; i<height; i++)
	{
		in_image_data[i] = new double*[width];
		for (j = 0; j<width; j++)
			in_image_data[i][j] = new double[3];
	}

	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			in_image_data[i][j][0] = luv_data[(i*width * 3) + j * 3];	//L
			in_image_data[i][j][1] = luv_data[(i*width * 3) + j * 3 + 1];	//u
			in_image_data[i][j][2] = luv_data[(i*width * 3) + j * 3 + 2];	//v
		}

	// Scale Luv values into [0,255]
	scaled_luv_data = scale(luv_data);

	// image containing the L component only
	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			l_data[(i*width * 3) + j * 3] = scaled_luv_data[(i*width * 3) + j * 3];
			l_data[(i*width * 3) + j * 3 + 1] = scaled_luv_data[(i*width * 3) + j * 3];
			l_data[(i*width * 3) + j * 3 + 2] = scaled_luv_data[(i*width * 3) + j * 3];
		}
	in_L_image = new wxImage(width, height, l_data);

	// image containing the u component only
	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			u_data[(i*width * 3) + j * 3] = scaled_luv_data[(i*width * 3) + j * 3 + 1];
			u_data[(i*width * 3) + j * 3 + 1] = scaled_luv_data[(i*width * 3) + j * 3 + 1];
			u_data[(i*width * 3) + j * 3 + 2] = scaled_luv_data[(i*width * 3) + j * 3 + 1];
		}
	in_u_image = new wxImage(width, height, u_data);

	// image containing the v component only
	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			v_data[(i*width * 3) + j * 3] = scaled_luv_data[(i*width * 3) + j * 3 + 2];
			v_data[(i*width * 3) + j * 3 + 1] = scaled_luv_data[(i*width * 3) + j * 3 + 2];
			v_data[(i*width * 3) + j * 3 + 2] = scaled_luv_data[(i*width * 3) + j * 3 + 2];
		}
	in_v_image = new wxImage(width, height, v_data);

}

unsigned char *ImageOperations::scale(double *luv_vector)
{
	int i, j, k;
	unsigned char *t;
	double max[3] = { luv_vector[0], luv_vector[1], luv_vector[2] };
	double min[3] = { luv_vector[0], luv_vector[1], luv_vector[2] };

	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
			for (k = 0; k<3; k++)
			{
				if (luv_vector[(i*width * 3) + j * 3 + k] < min[k])
					min[k] = luv_vector[(i*width * 3) + j * 3 + k];
				else if (luv_vector[(i*width * 3) + j * 3 + k] > max[k])
					max[k] = luv_vector[(i*width * 3) + j * 3 + k];
			}

	t = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));
	for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
			for (k = 0; k<3; k++)
			{
				t[(i*width * 3) + j * 3 + k] =
					(unsigned char)((luv_vector[(i*width * 3) + j * 3 + k] - min[k])
					* (min[k] != max[k] ? 255 / (max[k] - min[k]) : 0));
			}
	return t;
}
/* convert from CIE-L*u*v* colorspace to RGB colorspace
*/
double *ImageOperations::LuvToRGB(double *luv_pixel)
{
	double *rgb_pixel;
	double *xyz_pixel;
	double u0, v0;
	double uV, vV;	// u', v'
	rgb_pixel = new double[3];
	xyz_pixel = new double[3];

	// CIE-L*u*v* -> CIE-XYZ
	// Compute u0, v0 (corresponding to white color)
	u0 = 4 * 242.36628 / (242.36628 + 15 * 254.999745 + 3 * 277.63227);
	v0 = 9 * 254.999745 / (242.36628 + 15 * 254.999745 + 3 * 277.63227);

	uV = luv_pixel[1] / (13 * luv_pixel[0]) + u0;
	vV = luv_pixel[2] / (13 * luv_pixel[0]) + v0;

	// Y component
	xyz_pixel[1] = (pow(((double)(luv_pixel[0] + 16.0) / 116.0), 3.0))*254.999745;

	// X component
	xyz_pixel[0] = (-9 * xyz_pixel[1] * uV) / ((uV - 4)*vV - uV*vV);

	// Z component
	xyz_pixel[2] = (9 * xyz_pixel[1] - 15 * vV*xyz_pixel[1] - vV * xyz_pixel[0]) / (3.0*vV);

	// CIE-XYZ to RGB

	// R component
	rgb_pixel[0] = (xyz_pixel[0] * 3.240479 + xyz_pixel[1] * -1.537150 + xyz_pixel[2] * -0.498535);
	// G component
	rgb_pixel[1] = (xyz_pixel[0] * -0.969256 + xyz_pixel[1] * 1.875992 + xyz_pixel[2] * 0.041556);
	// B component
	rgb_pixel[2] = (xyz_pixel[0] * 0.055648 + xyz_pixel[1] * -0.204043 + xyz_pixel[2] * 1.057311);

	return rgb_pixel;
}


/* Create and display the output image based on the current labeling.
* Executed at each iteration.
*/
void ImageOperations::CreateOutput()
{
	int i, j;
	//unsigned = "0-255"
	unsigned char *out_data;

	/* Do not count GUI overhead
	*/
	timer.Stop();

	//allocates (width * height * 3 [possibly for RGB values]) memory blocks the size of an unsigned char 
	out_data = (unsigned char *)malloc(width*height * 3 * sizeof(unsigned char));

	// create output image
	for (i = 0; i<height; ++i)
		for (j = 0; j<width; ++j)
		{
			//stores pixel R channel in respective out_data slot
			out_data[(i*width * 3) + j * 3] =
				(unsigned char)out_regions[classes[i][j] * 3];
			//store pixel G channel value in respective out_data slot
			out_data[(i*width * 3) + j * 3 + 1] =
				(unsigned char)out_regions[classes[i][j] * 3 + 1];
			//stores pixel B channel in respective out_data slot
			out_data[(i*width * 3) + j * 3 + 2] =
				(unsigned char)out_regions[classes[i][j] * 3 + 2];
		}

	//previously frees allocated memory
	free(out_image);

	//creates a new image w/ out_data values
	out_image = new wxImage(width, height, out_data);

	// and display it
	((MyFrame *)frame)->GetOutputWindow()->SetScrollbars(10, 10, (out_image->GetWidth()) / 10, (out_image->GetHeight()) / 10);
	((MyFrame *)frame)->GetOutputWindow()->SetBmp(out_image);
	((MyFrame *)frame)->GetOutputWindow()->Refresh();
	((MyFrame *)frame)->GetOutputWindow()->Update();
	frame->RefreshRect(wxRect(645, 360, 100, 100));
	frame->Update();
	/* Continue timer
	*/
	timer.Start();
}


/* Metropolis & MMD

The third method (MMD) is a modified version of the Metropolis algorithm: at each iteration the new state is chosen randomly, but the decision to accept it is purely deterministic.
This is also a suboptimal technique but it is much faster than stochastic relaxation.
*/
void ImageOperations::Metropolis(bool mmd)
{
	//maps image pixels
	InitOutImage();

	CreateOutput();
	//wxString image_name;
	//image_name = "C:/Users/Elisabeth/Desktop/research-project/ColorMRFdemo/images";
	//SaveBmp(image_name);
	int i, j;
	int r;
	double kszi = log(alpha); // This is for MMD. When executing
	// Metropolis, kszi will be randomly generated.
	double summa_deltaE;

	TRandomMersenne rg(time(0));  // create instance of random number generator

	K = 0;
	T = T0; //set in the application window
	E_old = CalculateEnergy();

	do
	{
		summa_deltaE = 0.0;

		//iterating through each pixel
		for (i = 0; i<height; ++i)
			for (j = 0; j<width; ++j)
			{
				/* Generate a new label different from the current one with
				* uniform distribution.
				*/
				//if there are only two regions, set the curent region to 1-the region @ the current pixel
				if (no_regions == 2)
					r = 1 - classes[i][j];
				else //set the current region to a random region
					r = (classes[i][j] +
					(int)(rg.Random()*(no_regions - 1)) + 1) % no_regions;

				//checking to see if MMD should be run or not
				if (!mmd)  // Metropolis: kszi is a  uniform random number -- is it random??
					kszi = log(rg.Random());

				/* Accept the new label according to Metropolis dynamics.
				*/
				//if the value os kszi is less than or equal to local energy @ the current pixel minus 
				//the local energy at the pixel in the previously generated random region divided by the temperature...
				if (kszi <= (LocalEnergy(i, j, classes[i][j]) -
					LocalEnergy(i, j, r)) / T)
				{
					//add the absolute value of the local energy @ pixel w/ random region label - local energy @ pixel w/ predefined region label
					summa_deltaE +=
						fabs(LocalEnergy(i, j, r) - LocalEnergy(i, j, classes[i][j]));

					//E and e_old are both now equal to the old energy minus (local energy @ pixel w/ predefined region label +  local energy @ pixel w/ random region label)
					E_old = E = E_old -
						LocalEnergy(i, j, classes[i][j]) + LocalEnergy(i, j, r);

					//set the pixel's old region to the new region
					classes[i][j] = r;
				}
			}

		T *= c;         // decrease temperature
		++K;	      // advance iteration counter
		CreateOutput(); // display current labeling
	} while (summa_deltaE > t); // stop when energy change is small
}


/* ICM
*/
void ImageOperations::ICM()
{
	InitOutImage();
	CreateOutput();
	wxString image_name;
	image_name = "C:/Users/Ve/Documents/MS/thesis/ColorMRFdemo/images/output.bmp";
	SaveBmp(image_name);
	int i, j;
	int r;
	double summa_deltaE;

	K = 0;
	E_old = CalculateEnergy();

	do
	{
		summa_deltaE = 0.0;
		for (i = 0; i<height; ++i)
			for (j = 0; j<width; ++j)
			{
				for (r = 0; r<no_regions; ++r)
				{
					if (LocalEnergy(i, j, classes[i][j]) > LocalEnergy(i, j, r))
					{
						classes[i][j] = r;
					}
				}
			}
		E = CalculateEnergy();
		summa_deltaE += fabs(E_old - E);
		E_old = E;

		++K;	      // advance iteration counter
		CreateOutput(); // display current labeling
	} while (summa_deltaE > t); // stop when energy change is small
}


/* Gibbs sampler
*/
void ImageOperations::Gibbs()
{
	InitOutImage();
	CreateOutput();
	wxString image_name;
	image_name = "C:/Users/Ve/Documents/MS/thesis/ColorMRFdemo/images/output.bmp";
	SaveBmp(image_name);
	int i, j;
	double *Ek;		       // array to store local energies
	int s;
	double summa_deltaE;
	double sumE;
	double z;
	double r;

	TRandomMersenne rg(time(0)); // make instance of random number generator

	Ek = new double[no_regions];

	K = 0;
	T = T0;
	E_old = CalculateEnergy();

	do
	{
		summa_deltaE = 0.0;
		for (i = 0; i<height; ++i)
			for (j = 0; j<width; ++j)
			{
				sumE = 0.0;
				for (s = 0; s<no_regions; ++s)
				{
					Ek[s] = exp(-LocalEnergy(i, j, s) / T);
					sumE += Ek[s];
				}
				r = rg.Random();	// r is a uniform random number
				z = 0.0;
				for (s = 0; s<no_regions; ++s)
				{
					z += Ek[s] / sumE;
					if (z > r) // choose new label with probabilty exp(-U/T).
					{
						classes[i][j] = s;
						break;
					}
				}
			}
		E = CalculateEnergy();
		summa_deltaE += fabs(E_old - E);
		E_old = E;

		T *= c;         // decrease temperature
		++K;	      // advance iteration counter
		CreateOutput(); // display current labeling
	} while (summa_deltaE > t); // stop when energy change is small

	delete Ek;
}

//This function executes the Bee Algorithm
void ImageOperations::BeesAlg() {
	InitOutImage();

	//you can define a presegmentation output and send it here
	//CreateOutput();
	//wxString image_name;
	//image_name = "C:/Users/Ve/Documents/MS/thesis/ColorMRFdemo/images/output.bmp";
	//SaveBmp(image_name);
	//K = 0;
	E = 0;
	//number of scouts set to 2% of number of pixels
	int numberOfScout = height*width*0.02;
	//best sites set to 20% of number of scouts
	int bestSites = numberOfScout*0.2;
	//site size to be searched around best sites
	int siteSize = 2;
	scoutWithPollen = new Flower[numberOfScout];
	T = T0;

	TRandomMersenne rg(time(0));

	bool firstRun = true;

	//old energy
	E_old = CalculateEnergy();
	double summa_deltaE;
	do
	{

		summa_deltaE = 0.0;

		//if first time set all scouts randomly, else keep the best sites and the rest place randomly
		if (firstRun) {
			//there are no best sites on the first run
			placeRandomlyScouts(numberOfScout, 0);
			firstRun = false;
		}
		else {
			placeRandomlyScouts(numberOfScout, bestSites);
		}

		//ealuate pollen
		evaluateFitness(numberOfScout);

		//sort scouts based on pollen/fitness/less energy
		sortScouts(numberOfScout);
		//do neighborhood search
		siteSearch(numberOfScout, bestSites, siteSize);

		//new energy after search done
		E = CalculateEnergy();
		summa_deltaE += fabs(E_old - E);
		E_old = E;

		//T *= c;         // decrease temperature
		++K;	      // advance iteration counter
		CreateOutput(); // display current labeling

	} while (summa_deltaE > t);
}

void ImageOperations::placeRandomlyScouts(int n, int start) {
	TRandomMersenne rg(time(0));
	//n is the number of scouts
	int count = start;

	//continue the loop until the remaining scouts are no more
	while (count < n) {
		//generates random numbers
		double rh = rg.Random();
		double rw = rg.Random();
		int irh = rh*height;
		int irw = rw*width;

		Flower flower;
		flower.x = irh; //sets pixel for each flower
		flower.y = irw;
		bool exists = scoutExists(flower, count); //check to see if the current flower is available

		//if there isn't a scout at the current flower
		if (!exists) {
			flower.pollen = 0;
			scoutWithPollen[count] = flower;
			count++;
		}
	}
}

//checks to see if there is a scount at the randomly selected pixel chosen in placeRandomlyScouts
bool ImageOperations::scoutExists(Flower flower, int n) {
	bool found = false;

	for (int i = 0; i < n; i++) {
		Flower currentFlower = scoutWithPollen[i];
		if (currentFlower.x == flower.x && currentFlower.y == flower.y) {
			found = true;
			break;
		}
	}
	return found;
}

void ImageOperations::evaluateFitness(int n) {

	for (int i = 0; i < n; i++) {

		for (int r = 0; r < no_regions; r++) {
			//calculating the current energy at the scout's placement
			double existingLocalEnergy = LocalEnergy(scoutWithPollen[i].x, scoutWithPollen[i].y, classes[scoutWithPollen[i].x][scoutWithPollen[i].y]);

			//calculating the energy IF the scout were in the same location with a different label
			double localEnergy = LocalEnergy(scoutWithPollen[i].x, scoutWithPollen[i].y, r);

			//if a different label gives better energy, change the current label
			if (localEnergy < existingLocalEnergy) {
				classes[scoutWithPollen[i].x][scoutWithPollen[i].y] = r;
				scoutWithPollen[i].pollen = localEnergy;
			}
		}

	}
}

//clear pollen amount from scoutwithpollen array and return how many pixels's class is updated
void ImageOperations::sortScouts(int n) {

	for (int i = 0; i < (n - 1); i++)
	{
		for (int j = 0; j < n - i - 1; j++)
		{
			if (scoutWithPollen[j].pollen < scoutWithPollen[j + 1].pollen)
			{

				Flower temp = scoutWithPollen[j];
				scoutWithPollen[j] = scoutWithPollen[j + 1];
				scoutWithPollen[j + 1] = temp;
			}
		}
	}
}


void ImageOperations::siteSearch(int n, int bestSites, int siteSize) {

	for (int i = 0; i < bestSites; i++) {
		int x = scoutWithPollen[i].x;
		int y = scoutWithPollen[i].y;

		int startX, endX, startY, endY;

		if ((x - siteSize) < 0) {
			startX = 0;
		}
		else {
			startX = x - siteSize;
		}

		if ((x + 1 + siteSize) > height) {
			endX = height;
		}
		else {
			endX = x + 1 + siteSize;
		}

		if ((y - siteSize) < 0) {
			startY = 0;
		}
		else {
			startY = y - siteSize;
		}

		if ((y + 1 + siteSize) > width) {
			endY = width;
		}
		else {
			endY = y + 1 + siteSize;
		}


		for (int xi = startX; xi < endX; xi++) {

			for (int yi = startY; yi < endY; yi++) {

				for (int r = 0; r < no_regions; r++) {
					double existingLocalEnergy = LocalEnergy(xi, yi, classes[xi][yi]);

					double localEnergy = LocalEnergy(xi, yi, r);
					if (localEnergy < existingLocalEnergy) {
						classes[xi][yi] = r;


						scoutWithPollen[i].pollen = 0;
					}
				}

			}
		}
	}
}

/*

create initial population
keep track of ALL chromosomes -- old and new in the chromosome variable
the first [0,x] will be the unmutated chromosomes
the next [x,2*x] will be the mutated chromosomes

crossover and mutate the initial population -- only [0,x]
evaulate fitness of the entire chromosome variable (including mutated chromosomes)
pick solution
pick 25 chromosomes from [0,x] and 25 chromosomes from [0,2*x] for the new population
current challenge -- how to clear/reorder the 2d array?

address problem with last iteration

*/

void ImageOperations::GeneticAlgorithm() {
	InitOutImage();
	//debug
	//std::cout << CalculateEnergy() << std::endl;

	makeChromosomes();

	int num_it = 0; //iteration counter
	int max_it = numberOfIterations; //max number of iterations
	int siteSize = 5; //size of square radius
	bool lastIt = false; //checks for last iteration

	solution = new int*[height]; //stores solution
	std::copy(classes, classes + height, solution); //intial solution is the current segemntation

	do
	{
		if (num_it > 0) { //(num_it + 1) == max_it) {
			lastIt = true;
		}

		evaluate_fitness(siteSize, lastIt); //sends square radius and whether this is last iteration
		crossover();
		mutation();

		num_it++;
		K++;

		//uncomment this to use an ordered population in future iterations
		//formNewPopulation();

		std::copy(solution, solution + height, classes); //save the current solution
		chromosome_energy.clear();
		CreateOutput(); //display current labeling
	} while (num_it < max_it);

	//delete solution;
}

void ImageOperations::makeChromosomes() {

	int i, j, rndx, rndy, it;
	//chromosomes = new Chromosome*[no_chromosomes];
	chromosomes = new Chromosome*[no_chromosomes * 2]; //a 1D array of chromosomes -- height is equal to number of chromosomes + accounts for mutated chromosomes

	for (i = 0; i < (no_chromosomes * 2); i++) {
		chromosomes[i] = new Chromosome[no_regions]; //reserves n columns in each already-defined row where n = number of classes

		for (j = 0; j < no_regions; j++) { //in each spot in the array, put in a blank chromosome
			Chromosome null;
			null.val = -1;
			chromosomes[i][j] = null;
		}
	}

	//only fills in potential chromsome clusters for clearly-defined chromosomes -- NOT mutated chromosomses
	for (i = 0; i < no_chromosomes; i++) { //continue this loop for each row
		it = no_regions;
		do {


			rndx = rand() % height;
			rndy = rand() % width;

			int curr_label = classes[rndx][rndy]; //grabs random pixel label

			if (chromosomes[i][curr_label].val == -1) { //if this row's column is empty, then define this chromosome
				Chromosome temp;
				temp.x = rndx;
				temp.y = rndy;
				temp.val = classes[rndx][rndy];

				chromosomes[i][curr_label - 1] = temp;
				it--; //one less class to define in this current row of chromosomes
			}
		} while (it != 0); //continue this process until all the columns of this current row is filled
	}
}

void ImageOperations::formNewPopulation() {

	int i, j, curr_min_index;
	double curr_min;

	//use chromsome energy to find lowest energy points and their index in the array
	reordered_chromosomes = new Chromosome*[no_chromosomes * 2];

	for (i = 0; i < (no_chromosomes * 2); i++) {
		reordered_chromosomes[i] = new Chromosome[no_regions]; //reserves n columns in each already-defined row where n = number of classes

		for (j = 0; j < no_regions; j++) { //in each spot in the array, put in a blank chromosome
			Chromosome null;
			null.val = -1;
			reordered_chromosomes[i][j] = null;
		}
	}


	//find the smallest energy in chromosome energy
	for (i = 0; i < no_chromosomes; i++) {

		curr_min = chromosome_energy[0];
		curr_min_index = 0;

		for (int j = 0; j < chromosome_energy.size(); j++) {
			if (chromosome_energy[j] < curr_min)  {
				curr_min = chromosome_energy[j];
				curr_min_index = j;
			}
		}

		reordered_chromosomes[i] = chromosomes[curr_min_index];
	}

	chromosomes = reordered_chromosomes;

}

void ImageOperations::evaluate_fitness(int siteSize, bool firstIt) {

	int** classes_clone = new int*[height];

	std::copy(classes, classes + height, classes_clone); //create a copy of the current segmented image

	for (int i = 0; i < no_chromosomes * 2; i++) { //for each row of chromosomes

		std::copy(classes_clone, classes_clone + height, classes); //start with the original segmented image

		for (int j = 0; j < no_regions; j++) { //for each column in the current row

			int x = chromosomes[i][j].x;
			int y = chromosomes[i][j].y;

			int startX, endX, startY, endY;

			if ((x - siteSize) < 0) { //creating a square around the cluster center
				startX = 0;
			}
			else {
				startX = x - siteSize;
			}

			if ((x + 1 + siteSize) > height) {
				endX = height;
			}
			else {
				endX = x + 1 + siteSize;
			}

			if ((y - siteSize) < 0) {
				startY = 0;
			}
			else {
				startY = y - siteSize;
			}

			if ((y + 1 + siteSize) > width) {
				endY = width;
			}
			else {
				endY = y + 1 + siteSize;
			}

			for (int xi = startX; xi < endX; xi++) { //in each column and row surrounding the cluster center...

				for (int yi = startY; yi < endY; yi++) {

					//find which region provides the minimum energy and store that value in the pixel
					for (int r = 0; r < no_regions; r++) {
						double existingLocalEnergy = LocalEnergy(xi, yi, classes[xi][yi]);
						double localEnergy = LocalEnergy(xi, yi, r);

						if (localEnergy < existingLocalEnergy) {

							classes[xi][yi] = r;

						}
					}
				}
			}
		}

		//calculate energy for the newly segmented image
		chromosome_energy.push_back(CalculateEnergy());

		if (firstIt) { //i == 0) {
			//if this is the first iteration, the max and min values will be equal to the first energy calculation
			E_max = E_min = CalculateEnergy();
		}

		if (/*(i != 0)*/ !firstIt && CalculateEnergy() > E_max) {
			//if this is not the first iteration and the current energy is greater than the max energy, store this value
			E_max = CalculateEnergy();
		}

		if (!firstIt && CalculateEnergy() < E_min) {
			//if the current energy calculation is less than the min value, store this value
			E_min = CalculateEnergy();

			//store this segmented image as the solution as well
			std::copy(classes, classes + height, solution);
		}
	}

	//put the original segmented image back in the classes variable 
	std::copy(classes_clone, classes_clone + height, classes);
	//delete classes_clone;
}

void ImageOperations::crossover() {
	int half = no_regions / 2; //if width is odd, will round down due to int division
	int i, j;
	double temp;

	//single point crossover
	for (i = 0; i < no_chromosomes; i++) { //for each row
		for (j = half; j < no_regions; j++) {
			if ((chromosome_energy[i] / E_max) > c_prob) { //if the energy / maximum energy > crossover probability
				if (i == (no_chromosomes - 1)) {
					std::swap(chromosomes[i][j], chromosomes[0][j]);
					temp = chromosome_energy[i];
					chromosome_energy[i] = chromosome_energy[0];
					chromosome_energy[0] = temp;
				}
				else {
					std::swap(chromosomes[i][j], chromosomes[i + 1][j]);
					temp = chromosome_energy[i];
					chromosome_energy[i] = chromosome_energy[i + 1];
					chromosome_energy[i + 1] = temp;
				}
			}
		}
	}
}

/*
if there are 50 chromosomes,
chromosomes variable has 100 spaces available
0-49 = unmutated chromosomes
50-99 = mutated chromosomes
*/

void ImageOperations::mutation() {

	int i, j, rndx, rndy;
	int rnd_class, curr_label;

	/*
	copy the current chromosomes
	mutate the batch
	evaluate fitness -- choose 25
	*/

	//overload the copy function?
	//std::copy(chromosomes, no_chromosomes + height, mutated_chromosomes); //create a copy of the generated chromosomes

	//NO MORE MUTATION PROBABILITY!!!!
	for (i = 0; i < no_chromosomes; i++) { //for each row
		//if ((chromosome_energy[i] / E_max) < m_prob) { //if the energy/maximum energy < mutation probability...

		//copies over the current unmutated chromosome row to the mutated chromosome section
		for (j = 0; j < no_regions; j++) {
			chromosomes[i + no_chromosomes][j].x = chromosomes[i][j].x;
			chromosomes[i + no_chromosomes][j].y = chromosomes[i][j].y;
		}

		rnd_class = rand() % no_regions; //randomly choose a class
		curr_label = -1;

		while (curr_label != rnd_class) { //do this while the rnd_class is not equal to the current label
			rndx = rand() % height;
			rndy = rand() % width;

			curr_label = classes[rndx][rndy]; //picks a random pixel label 
		}

		chromosomes[i + no_chromosomes][rnd_class].x = rndx; //sets a new cluster center location
		chromosomes[i + no_chromosomes][rnd_class].y = rndy;

		//}
	}
}

