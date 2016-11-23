//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::PMTwf_ana+;
#pragma link C++ class PEhit+;
#pragma link C++ class std::vector<PEhit>+;
#pragma link C++ class std::vector< std::vector<PEhit> >+;
#pragma link C++ class wfInfo+;
#pragma link C++ class std::vector<wfInfo>+;
//ADD_NEW_CLASS ... do not change this line
#endif

