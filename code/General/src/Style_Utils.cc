#include "TH1.h"
#include "TPaveStats.h"

#include "Style_Utils.h"

void Move_StatBox(TH1* histo){
	TPaveStats *st = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
	st->SetX1NDC(0.65); //new x start position
	st->SetX2NDC(0.85); //new x end position
	st->SetY1NDC(0.6); //new x start position
	st->SetY2NDC(0.9); //new x end position
}				
