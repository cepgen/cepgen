#include "EventWriter.h"

EventWriter::EventWriter(const OutputType& type, const char* filename) :
  fType(type),
  fHepMCHandler(0)
{
  switch (fType) {
    case HepMC: { fHepMCHandler = new OutputHandler::HepMCHandler(filename); } break;
    default: return;
  }
}

EventWriter::~EventWriter()
{
  // HepMC persistent objects
  if (fHepMCHandler) delete fHepMCHandler;
}

void
EventWriter::operator<<(const Event* evt)
{
  switch (fType) {
    case HepMC: { (*fHepMCHandler) << evt; } break;
    default: return;
  }
}

