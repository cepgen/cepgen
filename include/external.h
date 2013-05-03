extern "C" {
  void grndm_(int&);
  void lujoin_(int&,int&); // connects a number of previously defined partons into a string configuration
  void luexec_();     // administrates the fragmentation and decay chain
  void lulist_(int&); // lists an event, jet or particle data, or current parameter values
  void luhepc_(int&); // converts between the LUJETS event record and the HEPEVT event record
}
