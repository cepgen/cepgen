#ifndef _EVENTSLIST_H
#define _EVENTSLIST_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "event.h"
#include "lheutils.h"

/**
 * Class containing all the Event objects generated in this run
 * @brief Collection of events in the run
 */
class EventsList {
 public:
  //EventsList(std::ofstream* of_, const HEPRUP *heprup_, const int dumpevery_=1000);
  EventsList(std::ofstream* of_, const int dumpevery_=1000);
  ~EventsList();
  //void AddEvent(Event *ev_, float weight_=1.) { this->_ev.push_back(std::pair<Event,float>(*ev_,weight_)); };
  //EventsList& operator=(const EventsList&);
  void AddEvent(Event *ev_);
  void StoreLHE(std::ofstream *of_);
  void Info();
  int NumEvents() { return this->_ev.size(); };
 private:
  void DumpEvents();
  //std::vector<std::pair<Event,float> > _ev;
  std::vector<Event> _ev;
  std::ofstream *_lheof;
  std::stringstream _lhe;
  int _dump_every;
};


#endif
