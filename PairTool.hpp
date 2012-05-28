#ifndef PAIRTOOL_H
#define PAIRTOOL_H

//Convert coorditantes into something.

class PairTool : public Main
{

protected:

  //Common variables between multiple frames.

  // dictionary for molecular name and properties
  map<string, FlexibleHandle> moldict;

  // Standard unit <=> internal unit
  Unit unit;

  // Interaction truncation
  Truncation rc;

  Box* box;

public:
  PairTool() : Main(), unit(), box(0)
  {
    cout.precision(17);
  }
  MolCollection* Read( FILE* );

  virtual ~PairTool()
  {
    if ( box != 0 ){
      delete box;
    }
    cerr << "~PairTool" << endl;
  }
};

#endif
